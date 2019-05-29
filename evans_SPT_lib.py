# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 15:53:29 2019

@author: Nate Evans

This is a supplementary package for methods developed in the Galbraith Lab [BME OHSU 2019] and is meant to be used with Trackpy [https://soft-matter.github.io/trackpy/v0.3.2/] and/or u-track.

"""

from skimage import io
import numpy as np
import scipy.misc
from matplotlib import pyplot as plt
import seaborn as sbn
from scipy.ndimage.filters import gaussian_filter
from skimage.morphology import disk
from skimage.filters import threshold_otsu, rank
from mpl_toolkits.mplot3d import Axes3D
import time
import pims
import pandas as pd
from scipy.io import loadmat
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
import statsmodels.formula.api as smf
from sklearn.cluster import DBSCAN
import trackpy as tp

class tif_img:
    '''
    This represents a greyscale tif stack acquired by TIRF imaging
    '''

    def __init__(self, dir_path):

        self.imgSeq = pims.TiffStack(dir_path) # this is what Trackpy uses
        self.img = io.imread(dir_path, plugin='tifffile')
        # self.__normalize_img()
        self.shape = self.img.shape
        meta = dir_path.split('/')[-1].split('_')
        self.date = meta[0]
        self.coverslip = meta[1]
        self.exp = meta[2]
        self.name = meta[3]
        self.exposure = meta[4]
        self.frames = meta[5][0:-4]
        self.fname = dir_path.split('/')[-1][0:-4]

    def __normalize(self, img):
        '''
        normalize the image, scale pixels to between 0,1

        inputs
            img <np.array> img to be normalized

        outputs
            normalized image
        '''
        return np.interp(img, (img.min(), img.max()), (0,1))

    def crop(self, xy0, xy1):
        '''
        This method crops the image to the dimensions provided.

        input
            xy0 <tuple> : (x0, y0) the upper left corner of the image, measured in pixels
            xy1 <tuple> : (x1, y1) the lower right corner of the image, measured in pixels

        output
            None

        '''
        self.setZeroImgSeq(xy0, xy1)

        x0,y0 = xy0
        x1,y1 = xy1
        img2 = np.zeros( (self.img.shape[0], x1-x0, y1-y0), dtype='uint16' )
        for i in range(self.img.shape[0]):
            img2[i] = self.img[i][y0:y1,y0:y1]
        self.img = img2

    def setZeroImgSeq(self, xy0, xy1):
        '''
        This acts like a crop but doesn't change the size of the image, just sets the cropped area to zero.

        input
            xy0 <tuple> : (x0, y0) the upper left corner of the image, measured in pixels
            xy1 <tuple> : (x1, y1) the lower right corner of the image, measured in pixels

        output
            None
        '''
        x0,y0 = xy0
        x1,y1 = xy1

        for i in range(self.img.shape[0]):
            a = self.imgSeq[i][y0:y1,y0:y1]
            self.imgSeq[i][:,:] = 0
            self.imgSeq[i][y0:y1,y0:y1] = a

    def plot3D(self, i=0, downscale=0.3):
        '''
        This function generates a 3D plot representing the pixels in x-y and intensity as z. Uses cubic interpolation.

        inputs
            i <int> (default 0) : index in tif stack; image choice to plot
            downscale <float> (default 0.3) : percentage to downscale, should produce plot size of (1-downscale)*orig_dims

        outputs
            None
        '''
        # downscaling has a "smoothing" effect
        lena = scipy.misc.imresize(self.img[i], downscale, interp='cubic')

        # create the x and y coordinate arrays (here we just use pixel indices)
        xx, yy = np.mgrid[0:lena.shape[0], 0:lena.shape[1]]

        # create the figure
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.plot_surface(xx, yy, lena ,rstride=1, cstride=1, cmap=plt.cm.gray,
                linewidth=0)

        # show it
        plt.show()

    def __calc_SNR(self, img, gaus_sd = 1, local_rad = 15, offset=70, show_imgs=False):
        '''
        This function calculates the Signal to noise ration for a single image. Signal to noise is usually measured by: (i0 - ib) / sqrt(i0) where i0 is signal intensity, and ib is background intensity. To approximate this, we use otsu threshholding to separate background and signal peaks. We then use the mean value of the signal mask as i0 and the mean signal of the background mask as ib.

        inputs
            img <np.array> : image to prefore SNR analysis on
            gaus_sd <int> : standard deviation in the gaussian filter applied before otsu threshold
            local_rad <int> : otsu neighborhood
            offset <int> : pixel intensity offset, higher value will decrease the area of the signal mask.
            show_imgs <boolean> : if true, will plot the mask images, for use in diagnostics and parameter tuning

        outputs
            SNR <float> : Signal to noise ratio; (i0-ib) / sqrt(i0)

        '''
        img2 = gaussian_filter(img, gaus_sd)

        selem = disk(local_rad)

        local_otsu = rank.otsu(img2, selem)

        if(show_imgs):
            fig, axes = plt.subplots(2,2)
            axes = axes.flatten()
            axes[0].imshow(img)
            axes[1].imshow(img2)
            axes[2].imshow(img2 >= local_otsu + offset)
            plt.show()

        i0 = np.average(img2[img2 >= (local_otsu + offset)])
        ib = np.average(img2[img2 < (local_otsu + offset)])

        SNR = (i0 - ib) / np.sqrt(i0)

        return SNR

    def gen_SNR_hist(self, n=50, show_imgs=False):
        '''
        This function calculates the signal-to-noise ratio for all images in an instances tif stack and then plots a histogram of the resulting values.

        inputs
            n <int> : the number of images to use in the SNR histogram; The first n frames in the tif stack are used.
            show_imgs <boolean> : if true, plots the SNR masks. For use in parameter tuning and diagnostics, use low n in these cases.

         outputs
             None
        '''
        start = time.time()
        SNP = []
        for i,img in enumerate(self.img[0:n]):
            if(i%(max(int(n/10), 1)) == 0 and i != 0 ):
                print('%.1f%% complete, frames processed: %d, time elapsed: %.2fs' %((i / n)*100, i, time.time()-start))
            SNP.append(self.__calc_SNR(img, show_imgs=show_imgs))

        plt.hist(SNP)
        plt.show()

        plt.plot(np.array(SNP))
        plt.show()

# set JAVA_HOME and JAVA_JDK paths
#conda install -c conda-forge pyjnius
#pip install imglyb
#pip install pyimagej
def run_thunderSTORM(path):
    '''
    wrapper for the imageJ plugin thunderSTORM: SPT localization program. Currently in dev, use with care.

    inputs
        path <str> : not implemented yet.

    outputs
        None
    '''

    myJ = ij.init(r'\Users\Nate\Downloads\fiji-win64\Fiji.apppip')
    myJ.op().run("Run analysis", "filter=[Wavelet filter (B-Spline)] scale=2.0 order=3 detector=[Local maximum] connectivity=8-neighbourhood threshold=std(Wave.F1) estimator=[PSF: Integrated Gaussian] sigma=1.3 method=[Weighted Least squares] full_image_fitting=true fitradius=3 fixed_intensity=true expected_intensity=500:900 nmax=10 pvalue=1.0E-6 mfaenabled=true keep_same_intensity=false renderer=[No Renderer]")


def import_thunderSTORM_results(path, px_size=111):
    '''
    imports the results generated from the fiji plugin thunderSTORM.

    inputs
        path <str> : name of the file, expected to be stored in thunderSTORM_outputs (relative path)
        px_size <int> : size of a pixel (sq) in nanometers
    '''
    data = pd.read_csv('./thunderSTORM_outputs/%s' %(path))

    assert list(data.columns.values) == ['id','frame','x [nm]','y [nm]','sigma [nm]','intensity [photon]','offset [photon]','bkgstd [photon]','uncertainty_xy [nm]'], 'check thunderSTORM results column outputs, did you include id?'

    #data.columns.values[2:4] = ['x','y']
    data = data.assign(x=data['x [nm]']/px_size)
    data = data.assign(y=data['y [nm]']/px_size)

    return data

def calculate_track_features(tracks, seconds_per_frame = 14e-3, nm_per_pixel = 111, diagnostics=False):
    '''
    This function generates indiviudal track features including:
        - track_dist : sum of all track lengths; total displacement [nm]
        - track_pc (optional) : principle component analysis of tracks, quantify evection by PCs (track direction)
            - PC1 : major axis of movement [unitless; but technically px]
            - PC2 : minor axis of movement [unitless; but technically px]
            - lambda1 : major axis variance [variance; unitless]
            - lambda2 : minor axis variance [variance; unitless]
            - evection_mag = lambda1/lambda2 : quantify 'evective magnitude' of the track movement (ellipse track cloud) [unitless]
        - avg_speed (nm / s) : track_dist (nm) / track_length * frame_time (s) [nm/s]
        - max_speed (nm / s) : max(track_dist_n) [nm/s]
        - min_speed (nm / s) : min(track_dist_n) [nm/s]
        - diffusion coefficient [nm**2 / s]
    '''

    print('beginning feature extraction, are your mappings correct?')
    print('using pixel mapping of %d nm/px...' %nm_per_pixel)
    print('using time mapping of: %f sec/frame' % seconds_per_frame)

    unique_ids = list(set(tracks['particle']))
    id_len = len(unique_ids)
    features = pd.DataFrame(index = unique_ids, columns = ['length', 'track_dist', 'dir', 'PC1', 'PC2', 'lambda1', 'lambda2', 'evection_mag', 'avg_speed', 'max_speed', 'min_speed', 'diffusion_coef', 'D_aic', 'D_r2', 'D_p'])


    for i, id in enumerate(unique_ids):
        if (i % int(len(unique_ids)/5) == 0):
            print('calculating track features: %d / %d complete' %(i, id_len))
        track_dat = tracks.loc[tracks['particle'] == id]

        # set track length
        features.at[id, 'length'] = len(track_dat.index)

        # get total track displacement
        disp, max, min, avg, dir, D, aic, p, r2 = get_track_feats(track_dat, dtm=seconds_per_frame, pxm=nm_per_pixel, diagnostics=diagnostics)
        features.at[id, 'track_dist'] = disp*nm_per_pixel # nm
        features.at[id, 'avg_speed'] = avg*nm_per_pixel/seconds_per_frame # nm / second
        features.at[id, 'max_speed'] = max*nm_per_pixel/seconds_per_frame # nm / second
        features.at[id, 'min_speed'] = min*nm_per_pixel/seconds_per_frame # nm / second
        features.at[id, 'dir'] = [_dir*nm_per_pixel for _dir in dir] # (nm, nm) vector from first to last track pt
        features.at[id, 'diffusion_coef'] = D # nm^2 / sec
        features.at[id, 'D_aic'] = aic # aikike information criteria, unitless, smaller is better for comparing models
        features.at[id, 'D_p'] = p # p value calculated from the f-statisitc
        features.at[id, 'D_r2'] = r2 # R^2 value, % variability explained by the model.

        # get track percent aggregate
        features.at[id, 'percent_aggregate'] = np.sum(track_dat['is_agg']) / len(track_dat)

        # run PCA
        PC1, PC2, lambda1, lambda2 = get_track_pca(track_dat)
        features.at[id, 'PC1'] = PC1
        features.at[id, 'PC2'] = PC2
        features.at[id, 'lambda1'] = lambda1
        features.at[id, 'lambda2'] = lambda2
        features.at[id, 'evection_mag'] = lambda1 / lambda2

        features.Data_Dictionary = {'PCx' : 'x Principle Component of track [x, y; direction-unitless]; direction of major/minor axis of the track location data',
        'lambdaX' : 'eigenvalues for PCx; variance in PCx [unitless]',
        'evection_mag' : 'lambda1/lambda2; measure of evection magnitude or how elliptical the track location data is [unitless]',
        'track_dist': 'the total distance the track covers; sum of all displacements [nm]',
        'xxx_speed' : 'xxx track speed [nm / sec]',
        'dir' : 'track direction, as calculated by last_track - first_track, not normalized. [nm,nm]',
        'diffusion_coef' : 'Track Diffusion Coefficient; calculated as slope/4 of MSDs vs lagTime, using pmin ~ 1/4*N',
        'D_aic' : 'aikike information criteria for Diffusion coefficient fit [unitless]; a measure of model fit (smaller is better) but not comparable between models; for use comparing different model types (polynomial degree)',
        'D_p' : 'p-value of the diffuision coef fit as calculated from the F-statistic [probability of null hypothesis; probability model occured by chance]',
        'D_r2': 'R^2 Value for the diffusion coef fit; variability explained by the model [unitless]',
        'percent_aggregate' : 'the probability of each track being apart of an aggregate'}

    return features

def get_track_pca(dat):

    pca = PCA(n_components = 2)

    x = list(dat['x'])
    y = list(dat['y'])
    X = np.array([np.array([x_, y_]) for x_, y_ in zip(x,y)])

    pca.fit(X)

    eigenvects = pca.components_
    eigenvals = pca.explained_variance_

    #print(eigenvects)
    #print(eigenvals)

    return eigenvects[0], eigenvects[1], eigenvals[0], eigenvals[1]

def get_track_feats(dat, dtm, pxm, diagnostics=False):

    sum = 0
    speeds = []

    first = 0
    last = 0

    positions = []
    f_min = int(min(list(dat['frame'])))
    f_max = int(max(list(dat['frame'])))

    for i in range(f_min + 1, f_max + 1) :
        last = dat.loc[dat['frame'] == (i-1)]
        now =  dat.loc[dat['frame'] == i]

        lxy = ( last['x'].iloc[0], last['y'].iloc[0] )
        nxy = (now['x'].iloc[0], now['y'].iloc[0])

        # last frame to this frame
        dist = np.sqrt( (nxy[0]-lxy[0])**2 + (nxy[1] - lxy[1])**2 ) # px

        if (i == f_min + 1) :
            first = lxy
            positions.append(lxy)

        if (i == (f_max) ) :
            last = nxy

        positions.append(nxy)
        speeds.append(dist) # px / frame
        sum += dist # calc dist

    try:
        # alternative: https://github.com/soft-matter/trackpy/blob/master/trackpy/motion.py
        D, aic, p, r2 = calculate_diffusion_coefficients(positions, dat, dt=dtm, pxm=pxm, diagnostics=diagnostics)
    except:
        print('warning: tracks with less than 3 points can not calculate a diffusion coefficient; suggest filtering for stubs > 5')
        raise

    # direction (x1-x0, y1-y0)
    movement_direction = [last[0] - first[0], last[1] - first[1]]

    return sum, max(speeds), min(speeds), np.mean(speeds), movement_direction, D, aic, p, r2

def calculate_diffusion_coefficients(positions, dat, dt, pxm, diagnostics=False):
    '''
    We follow the method described in:
    [1]
        Michalet, Xavier. “Mean square displacement analysis of single-particle
        trajectories with localization error: Brownian motion in an isotropic
        medium.” Physical review. E, Statistical, nonlinear, and soft matter
        physics vol. 82,4 Pt 1 (2010): 041914. doi:10.1103/PhysRevE.82.041914

    and use suggestions found in section 3 of:
    [2]
        Michalet, Xavier, and Andrew J Berglund. “Optimal diffusion coefficient
        estimation in single-particle tracking.” Physical review. E,
        Statistical,] nonlinear, and soft matter physics vol. 85,6 Pt 1 (2012):
        061916. doi:10.1103/PhysRevE.85.061916

    Optimizations for pmin(x,N) are available in the appendix of [2] but for
    now, we will use the suggestions from below that pmin ~ 0.25*N
    [3]
        Saxton, M J. “Single-particle tracking: the distribution of diffusion
        coefficients.” Biophysical journal vol. 72,4 (1997): 1744-53.
        doi:10.1016/S0006-3495(97)78820-9

    MSD = pn = 1/(N-n) * sum(pos[i+n]-pos[i]) from i = 0:(N-n) where n = 1,2,3...,N-1

    p(t) = ρ(iΔt) = 2d(σ2 − 2RDΔt) + 2dDiΔt. [2]
        where
            d = dimensions of movement freedom
            D = Diffusion coefficient
            t = frame
            n = lag
            N = Total number of localizations

        by fitting a linear model to our measured p(t) we obtain:
            p(t) = a + bt
            where
                D = b/4
                a = 4(σ2 − 2RDΔt)

        for now, fit using pn[0:0.75*N]

        --- TODO:: pmin optimization ---
        "The main results of this work can be summarized as follows. In the presence of a localization error1 σ,
        the critical control parameter is the reduced localization error x = σ2/DΔt, where D is the diffusion
        constant and Δt is the frame duration.

        When this dimensionless ratio x ≪ 1, the best estimate of the diffusion coefficient is obtained using
        the first two points of the MSD curve (excluding the (0, 0) point).

        When x ≫ 1, the standard deviation of the first few MSD points is dominated by localization uncertainty,
        and therefore a larger number of MSD points are needed to obtain a reliable estimate of D. The optimal
        number pmin of MSD points to be used depends only on x and N, the number of points in the trajectory. For
        small N, the optimal number pmin of points may sometimes be as large as N, while for large N, pmin may be
        relatively small." [1]

        [2] offers an approximation of pmin(x, N) for D and sigma.

        [3] suggests that pmin ~ 1/4 * N <- this is what we will use for now

    inputs
        positions <list> track distances stored as a list of tuples: [(x1, y1,), (x2, y2)...(xn, yn)]. Units are in pixels.
        dt <float> frame time mapping, default is set at 14ms/frame.
        pxm <float> pixel mapping, default is set to 111nm/px
        diagnostics <boolean> whether to print diagnostic summary and plots

    outputs
        D <float> Diffusion Coefficient
        aic <int> aikike information criteria; measure of fit on MSD vs time OLS regression.
        p <float> probability of null hyp. on regression fit
        r2 <float> proportion of variance explained by regression fit
    '''

    pn = []
    lagTime = []
    N = len(positions)
    # n is the number of frames to measure between

    '''
    for n in range(1,N): # 1->(N-1) non-inclusive
        sum = 0
        for i in range(0, N-n):
            disp_px = ((positions[i+n][0] - positions[i][0])**2 + (positions[i+n][1] - positions[i][1])**2)**0.5     # euclidean distance
            disp_nm = disp_px*pxm
            sum += disp_nm**2
        pn.append(sum / (N-n)) # nm
        lagTime.append(dt*n) # seconds
    '''

    micronsperpixel = pxm/1000
    framespersecond = 1/dt
    #print('mpp: %.2f \nfps: %.2f' %(micronsperpixel, framespersecond))
    # dat is in pixels/frame, pn is microns/second
    pn = tp.msd(dat, mpp=micronsperpixel, fps=framespersecond, max_lagtime=100, detail=False, pos_columns=None)

    # pmin optimizations
    # see appendix of [2] for pmin approximation functions
    pmin = 1*N

    df = pd.DataFrame(columns=['y', 'x'])
    df['x'] = np.array(pn['lagt'])[0:int(pmin)]
    df['y'] = np.array(pn['msd'])[0:int(pmin)]

    weights_linear = np.polyfit(df['x'], df['y'], 1)
    weights_quadratic = np.polyfit(df['x'], df['y'], 2)

    linear_model = np.poly1d(weights_linear)
    quadratic_model = np.poly1d(weights_quadratic)

    res_lin = smf.ols(formula='y ~ x', data=df).fit()
    res_quad = smf.ols(formula='y ~ x + np.power(x,2)', data=df).fit()

    # QUADRATIC FIT : qa*x^2 + qb*x + qc
    qc = res_quad.params[0]
    qb = res_quad.params[1]
    qa = res_quad.params[2]
    qbic = res_quad.bic
    df = df.assign(quadfit = lambda v: qa*v.x**2 + qb*v.x + qc)

    # LINEAR FIT : lb*x + lc
    lb = res_lin.params[1]
    lc = res_lin.params[0]
    lbic = res_lin.bic
    df = df.assign(linfit = lambda v: lb*v.x + lc)

    if (lbic < qbic):
        if(diagnostics):
            print('LINEAR MODEL IS SMALLER BY: %.2f' %(lbic-qbic))
    else:
        if(diagnostics):
            print('QUADRATIC MODEL IS BEST BY: %.2f' %(qbic-lbic))

    '''
    model = sm.OLS(Y, X)
    results = model.fit()
    a = results.params[0] # nm^2
    b = results.params[1] # nm^2 / second
    D = b/4
    x = (a/4) / (D*dt)
    aic = results.aic
    p = results.f_pvalue
    r2 = results.rsquared
    '''

    if (diagnostics):
        #print(pn)
        #print('D, aic, p, r2, x, a, b: ' + str([D, aic, p, r2, x, a, b]))
        #print(res_lin.summary())
        #print(res_quad.summary())
        x1 = [xy[0] for xy in positions]
        y1 = [xy[1] for xy in positions]
        df_track = pd.DataFrame({'x': x1, 'y':y1, 'lab':range(len(x1))})
        fig, ax = plt.subplots(1,2, figsize=(15,10))
        sbn.lineplot(x='x', y='y', data =df_track, ax=ax[0], sort=False).set_title('Particle Track')
        [ax[0].text(x, y, lab, fontsize=7, color='red') for x,y,lab in zip(x1,y1,range(len(x1)))]

        sbn.scatterplot(x=pn['lagt'], y=pn['msd'], color='blue', ax=ax[1])
        sbn.lineplot(x='x', y='linfit', color='red', ax=ax[1], data=df)
        sbn.lineplot(x='x', y='quadfit', color='green', ax=ax[1], data=df).set_title('MSD vs LagTime')

        ax[0].set(xlabel='x [px]', ylabel='y [px]')
        ax[1].set(xlabel='lagtime [n*dt] seconds', ylabel='Mean Squared Displacement [um^2]')
        plt.show()

    return -1,-1,-1,-1#D, aic, p, r2

def detect_aggregations(tracks, n=15, _eps = 0.275, pxmap = 111, plot=False):
    '''
    Detects aggregation by points with at least n colocalizations. This is done by density based clustering of all tracks, without time separation.

    input
        tracks <dataframe> must have, at least: x,y,particle(id)
        n <int> number of localizations required in a cluster to be considered an aggregate
        plot <boolean> decision to plot cluster results

    output
        tracks <dataframe> tracks from input with new feature ['aggregate'], where 0 pertains to non-aggregate, and each unique int >0 corresponds to distinct aggregate.
    '''
    print('starting aggregate detection on %d localizations' %len(tracks.index))
    print('beginning clustering...')
    x = np.array(tracks['x'])
    y = np.array(tracks['y'])
    X = np.array([[_x,_y] for _x,_y in zip(x,y)])

    # cluster the bound particles

    print('Using eps of %.3f [%.2f nm] - does that sound right? ~30nm is default.' %(_eps, (_eps*pxmap)))
    agg = DBSCAN(eps = _eps, min_samples=n).fit_predict(X)

    tracks['aggregate'] = pd.Series(agg, index = tracks.index)

    tracks['is_agg'] = pd.Series( (tracks['aggregate'] != -1), index = tracks.index) # double check -1 is unlabeled

    if (plot):
        plt.figure(figsize=(15,15))
        sbn.scatterplot(x='x',y='y',hue='is_agg',alpha=0.05, data=tracks)
        plt.show()

    print('clustering complete. %d localizations labeled as aggregate.' % np.sum(tracks['is_agg']))
    return tracks



if __name__ ==  '__main__' :
    '''
    testing and diagnostics for development.
    '''

    print('beginning evans_SPT_lib __main__ ')

    TNT_OUTPUT_PATH = '/Users/galbraithlab/Desktop/evans/tutorial_tnt_outputs' #'./../tutorial_tnt_outputs'
    TNT_OUTPUT_NAME = '180920_Cs1C2_CAD_JFHaloLifeactRac_14ms_2Kfr_02-CROPPED_2019-m04-d07-18h59_TNT.mat'
        # ANALYSIS

    # load mat returns the data as a dictionary
    sptRes = loadmat(TNT_OUTPUT_PATH + '/' + TNT_OUTPUT_NAME)

    # id, frame, x, y, z, amplitude - for trackpy, let particle = id
    tracks = pd.DataFrame(sptRes['trackingData'], columns = ['particle', 'frame', 'x', 'y', 'z', 'amplitude'] )

    tracks2 = tp.filter_stubs(tracks, 10)
    tracks2 = tracks2.reset_index(drop=True)
    print(tracks.head())
    print(tracks2.head())

    tracks = detect_aggregations(tracks2, n=15, plot=True)

    # tracks must have AT LEAST 3 points;
    track_feats = calculate_track_features(tracks, seconds_per_frame = 14e-3, nm_per_pixel = 111, diagnostics=True)
    print(track_feats.head())

    # main
    #tif_path=r'\Users\Nate\Box Sync\Galbraith_Research\data\180920_Cs1C2_CAD_JFHaloLifeactRac_14ms_2Kfr_02-CROPPED.tif'

    #run_thunderSTORM(tif_path)

    #thunderSTORM_res_path = '180920_Cs1C2_CAD_JFHaloLifeactRac_14ms_2Kfr_02-CROPPED_thunderSTORM-localization.csv'
    #import_thunderSTORM_results(thunderSTORM_res_path)
