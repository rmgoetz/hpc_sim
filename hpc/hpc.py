# Mauricio Diaz-Ortiz 
# 27 March 2020
# mdiazort@ufl.edu

import pykat as pk
import numpy as np
import matplotlib.pyplot as plt

class CMOS_sensor:
    """
    This class will generate a CMOS object that will provide x/y arrays for generating images
    and then provide functionality for adding noise to that images as well as digitization.
    Default values are those for the Blackfly S, however this can be used generally.
    """
    def __init__(self, pixel_pitch=6.9e-6, x_resolution=720, y_resolution=540, exposure_time=4e-6, quantum_eff=0.03, pixel_well_depth=22187):
        "Intializes a CMOS sensor with requested dimensions and generates x and y arrays to generate desired images."
        self.pixel_pitch = pixel_pitch
        self.x_resolution = x_resolution
        self.y_resolution = y_resolution
        self.exposure_time = exposure_time
        self.quantum_eff = quantum_eff
        self.pixel_well_depth = pixel_well_depth
        self.x_array = np.linspace(-self.pixel_pitch*(self.x_resolution/2), self.pixel_pitch*(self.x_resolution/2), x_resolution)
        self.y_array = np.linspace(-self.pixel_pitch*(self.y_resolution/2), self.pixel_pitch*(self.y_resolution/2), y_resolution)

    def set_pixel_well_depth(pixel_well_depth):
        self.pixel_well_depth = pixel_well_depth

    def set_quantum_eff(quantum_eff): 
        self.quantum_eff = quantum_eff

    def set_exposure_time(exposure_time):
        self.exposure_time = exposure_time

    def add_read_noise(self, image, mean=3.71):
        "For a given mean will add noise representative of the electronics noise in the system"
        noise = np.abs(np.random.normal(loc=mean, size=image.shape))
        noisey_image = image + noise
        return noisey_image

    def add_shot_noise(self, image):
        shot_noise = np.random.poisson(np.real(image))
        return shot_noise

    def convert_to_photons(self, image, wavelength=1064e-9):
        h = 6.6260689158E-34 #J/s
        c = 2.99792458E8
        lam = 1064e-9
        E = h*c/lam #photon energy

        image_E = np.real(image)*self.pixel_pitch**2*self.exposure_time
        image_photons = np.floor(image_E/E)
        image_photons = image_photons.astype(int)
        return image_photons

    def convert_to_electrons(self, image, convert_to_photons=False):
        if convert_to_photons == False:
            image_electrons = image*self.quantum_eff
        else:
            image_photons = self.convert_to_photons(image)
            image_electrons = image_photons*self.quantum_eff

        return image_electrons

    def digitize(self, image, bitdepth):
    	"provided a given bitdepth and image will produce bins scaled to the max well depth of the sensor and place all vlaues in an appropriate bin. Must feed in an image that has been converted to electrons."
    	bits = int(2**bitdepth)
    	bins = np.linspace(0, self.pixel_well_depth, bits)
    	digitized_image = np.digitize(np.real(image), bins)
    	return digitized_image

    def capture(self, image, bitdepth, mean=3.71):
        "Given an intensity image, will produce a more realistic version as if passing through the camera."
        photons = self.convert_to_photons(image)
        shot_noise = self.add_shot_noise(photons)
        electrons = self.convert_to_electrons(shot_noise)
        read_noise = self.add_read_noise(electrons, mean)
        digitized_image = self.digitize(read_noise, bitdepth)
        return digitized_image

class beam:
    def __init__(self, power, w0, z, freq=None, spatial="gauss"):
        self.power = power
        self.w0 = w0
        self.freq = freq
        self.spatial = spatial
        self.z = z 

    def set_power(self, power):
        self.power = power

    def set_w0(self, w0):
        self.w0 = w0

    def set_z(self, z):
        self.z = z

    def set_freq(self, freq):
        self.freq = freq

    def set_spatial(self, spatial):
        if spatial == "gauss" or "flattop":
            self.spatial = spatial
        else:
            print('Entered spatial profile is not recognized. Please enter either gauss or flattop.')

    def generate_amplitude_map(self, x_array, y_array, pixel_pitch=0, array=None, max_val=None, x_offset=0, y_offset=0): #how will user input/upload array??
        "Given an x and y array will produce an amplitude map of the beam as defined. x/y offsets will displace the beam in the respective axis."
        if self.spatial == "gauss":
            import pykat.optics.gaussian_beams as gb
            q = gb.BeamParam(w0=self.w0, z=self.z)
            HG00 = gb.HG_mode(q,n=0, m=0)
            u00 = np.sqrt(self.power)*HG00.Unm(y_array-y_offset, x_array-x_offset)
            return u00
        
        elif self.spatial == "user":
            #applies user array to custimize beam shape
            yy = np.meshgrid(x_array,y_array)
            ampIntArr = np.zeros(np.shape(yy))#(ROWS, COLUMNS)
            
            if array.shape != (np.shape(yy)): #numpy array is (ROWS, COLUMNS)
                print('Entered array is not the correct size. Please enter an array with ' 
                    + str(np.shape(yy)) + ' Rows and Columns.')
            
            if array.dtype == bool:
                ampIntArr[array] = np.sqrt(self.power)/(pixel_pitch * np.sqrt(numPix)) 
                return ampIntArr
            else: #assuming the entered array isn't a string or anything stupid...
                if max_val == None:
                    max_val = array.max()
                normArr = np.zeros(np.shape(yy))
                numPix = np.count_nonzero(array) #number of pixels with nonzero value
                normArr = array/max_val
                ampIntArr = np.sqrt( normArr * self.power / (numPix * np.sum(normArr)) ) / pixel_pitch
                return ampIntArr

        else:
            xx, yy = np.meshgrid(x_array,y_array)
            tophat = np.zeros(np.shape(xx))
            r = 1e-3
            mask = xx**2 + yy**2 < r**2
            tophat[mask] = np.sqrt(self.power)/(np.sqrt(np.pi)*r)
            return tophat

def add_RIN(image, mean):
    "This will be made more complex in later versions, but will add random noise on top of the image"
    noise = np.abs(np.random.normal(loc=mean, size=image.shape))
    noisey_image = image + noise
    return noisey_image

def four_point(images):
    "Four point phase algorithm"
    phase = []

    if len(images) < 4:
        print("There are not enough images to calculate the phase!")

    else:
        for i in range(len(images)):
            if (i+1)%4 == 0 :
                I1 = images[i-3]
                I2 = images[i-2]
                I3 = images[i-1]
                I4 = images[i]

                phi = np.arctan2((I4-I2), (I1-I3))
                phase.append(phi)

    return phase

def carre(images):
    "Carre phase algorithm"

    phase = []
    if len(images) < 4:
        print("There are not enough images to calculate the phase!")

    else:
        for i in range(len(images)):
            if (i+1)%4 == 0 :
                I1 = images[i-3]
                I2 = images[i-2]
                I3 = images[i-1]
                I4 = images[i]
            #print(type(I1), type(I2),type(I3),type(I4))
                A = I2-I3
                B = I1-I4
                num = (A+B)*(3*A-B)
                num = np.sqrt(abs(num))
                pm = np.sign(A)
                #num = np.sign(I2-I3)*np.sqrt(np.abs(I1+I2-I3-I4)*(3*(I2-I3)-I1+I4))
                denom = I2+I3-I1-I4
                phi = np.arctan2(pm*num, denom)
                phase.append(phi)

    return phase  
    
def novak(images):
    'Novak phase algorithm'
    phase = []
    
    for i in range(len(images)):
        if i%4 == 0 and i != 0:
            I1 = images[i-4]
            I2 = images[i-3]
            I3 = images[i-2]
            I4 = images[i-1]
            I5 = images[i]
            
            d = 2*I3 - I1 - I5
            a24 = I2 - I4
            a15 = I1 - I5
            n = np.sqrt(abs(4*a24**2 - a15**2))
            pm = np.sign(a24)
            phase.append(np.arctan2(pm*n,d))

        else: continue
    
    return phase

def time_array(no_images, collection_time):
	'produces a time array resprestative of the spacing between images collected over a period of time'
	t = np.linspace(0, collection_time, no_images)
	return t

def generate_beatnote(amp_1, amp_2, bn_frequency, time_array):
	"Generates an array containing the expected beat note from two gaussian amplitudes"
    # add phase offset
	phase = [2*np.pi*bn_frequency*i for i in time_array]
	amp1_stepped = [amp_1*np.exp(1j*i) for i in phase]
	amp_total = [amp_2+i for i in amp1_stepped]
	I = amp_total*np.conj(amp_total)
	return I

def animate_images(images, cmap='gray', interval=25, cbar_lim=None):
    "Produces an animation given a set of images, interval is in ms."

    from matplotlib import animation
    from IPython.display import HTML

    fig = plt.gcf()
    im = plt.imshow(np.real(images[0]), cmap=cmap)
    fig.colorbar(im)

    if cbar_lim != None:
        plt.clim(vmin=0,vmax=cbar_lim)

    plt.close()
    
    def animate(i):
        im.set_data(np.real(images[i]))
        return im,
    
    anim = animation.FuncAnimation(fig, animate, frames=range(0,len(images)), interval=interval, repeat=True)
    display(HTML(anim.to_jshtml()))

def plot(image, cmap='gray', cbar_lim=None):
    plt.imshow(image, cmap=cmap)
    plt.colorbar()

    if cbar_lim != None:
        plt.clim(vmin=0, vmax=cbar_lim)

def gauss_amp(x_array, y_array, w0, z, x_offset, y_offset):
    import pykat.optics.gaussian_beams as gb
    q = gb.BeamParam(w0=w0, z=z)
    HG00 = gb.HG_mode(q,n=0, m=0)
    u00 = HG00.Unm(y_array-y_offset, x_array-x_offset)
    return u00
