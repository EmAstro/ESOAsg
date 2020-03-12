"""
Module to perform some simple operations on an image
"""





def findSources(imgData,
                imgStat,
                sigDetect=3.501,
                minArea=16.):
    """Detecting sources on image using sep:
    https://github.com/kbarbary/sep
    This is simply running SExtractor on the input
    image and returning relevant parameter. It returns
    location of each source in pixel, semi-major and
    semi-minor axis and rotation angle in degree. Only
    good sources are returned. By meaning sources that
    have SExtractor flag = 0. The full set of sources
    with all information is also returned in the
    astropy-table: allObjects.
    Note that no filtering is considered (i.e. filter_kernel=None)
    to avoid single pixels to be relevant into the
    detection.
    Parameters
    ----------
    imgData : np.array
        2D image of the field
    imgStat : np.array
        2D variance image from imgData
    sigDetect : np.float
        detection limits for source (default=3.5)
    minArea : np.float
        minimum area to consider a source
        as detected (default=16.)
    Returns
    -------
    xPix, yPix, aPix, bPix, angle : np.arrays
        position, dimension, and angle (in degrees)
        of the good sources
    allObjects : astropy table
        full info on detected sources
    """

    print("findSources: Starting sources detection")
    print("findSources: Creating background model")
    bgMedian = np.nanmedian(imgData)
    bgSigma = np.sqrt(np.nanmedian(imgStat))
    bgMask = np.zeros_like(imgData)
    bgMask[(np.abs(imgData-bgMedian)>7.*bgSigma)] = np.int(1)
    imgBg = sep.Background(imgData, mask=bgMask,
                           bw=64., bh=64., fw=5., fh=5.)
    imgDataNoBg = np.copy(imgData) - imgBg
    print("findSources: Searching sources {}-sigma above noise".format(sigDetect))
    allObjects = sep.extract(imgDataNoBg, sigDetect,
                             var=imgStat,
                             # var=np.nanmedian(imgStat),
                             minarea=minArea,
                             filter_type='matched',
                             gain=1.1,
                             clean=True,
                             deblend_cont=0.3,
                             filter_kernel=None)
    # Sorting sources by flux at the peak
    indexByFlux = np.argsort(allObjects['peak'])[::-1]
    allObjects = allObjects[indexByFlux]
    goodSources = allObjects['flag']<1
    xPix = np.array(allObjects['x'][goodSources])
    yPix = np.array(allObjects['y'][goodSources])
    aPix = np.array(allObjects['a'][goodSources])
    bPix = np.array(allObjects['b'][goodSources])
    angle = np.array(allObjects['theta'][goodSources])*180./np.pi
    print("findSources: {} good sources detected".format(np.size(xPix)))
    # Deleting temporary images to clear up memory
    del imgBg
    del imgDataNoBg
    del goodSources
return xPix, yPix, aPix, bPix, angle, allObjects
