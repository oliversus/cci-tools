import pandas as pd

def init():
    global data_folder, figuresDir, latex_variables, NA1, NA2, SIB, AFR, main_folder, sceneTime, sceneTimes, \
        sceneName, sceneNames, calipso_segments, boundaries_ENV, boundaries_MODIS
    main_folder = "/cmsaf/esa_doku/ESA_Cloud_cci/publications/CC4CL_paper/"
    data_folder = main_folder + "data/"
    figuresDir = main_folder + "figures/"
    latex_variables = dict()
    sceneTimes = ["07221915", "07222058", "07270810", "01071228", "10241345"]
    NA1 = sceneTimes[0]
    NA2 = sceneTimes[1]
    SIB = sceneTimes[2]
    TIM = sceneTimes[3]
    AFR = sceneTimes[4]
    sceneTime = ""
    sceneName = ""
    sceneNames = ['NA1', 'NA2', 'SIB', 'TIM', 'AFR']
    calipso_segments = pd.DataFrame([
        (60.0, 70, 75, 85.0), # NA1
        (50.0, 58.5, 66.0, 75.0), # NA2
        (55,71,76,88), # SIB
        (0,0,0,0), # TIM
        (-30,-7,12,30)], # AFR
        columns=list('1234'),
        index=sceneNames)
    boundaries_ENV = pd.DataFrame([
        (-129.,62,-63,83.5,-121,58.5,-68,78.5), # NA1
        (-133.5, 52., -115.5, 75., -129., 45., -105., 71), # NA2
        (36,62.5,110,83.5,45,60,105,78.5), # SIB
        (0,0,0,0,0,0,0,0), # TIM
        (0,0,0,0,0,0,0,0)], # AFR
        columns=list('12345678'),
        index=sceneNames)
    # values: eastern boundary (lon1, lat1, lon2, lat2), western boundary (lon1, lat1, lon2, lat2)
    boundaries_MODIS = pd.DataFrame([
        (-125,60,-177,72,-78,57,-70,80), # NA1
        (-132.5, 45., -168., 68, -102.5, 45., -100.5, 75.), # NA2
        (40,62.5,-10,72,91,58,97,80), # SIB
        (0,0,0,0,0,0,0,0), # TIM
        (-10,-13,-17,18,12.5,-15,6,18)], # AFR
        columns=list('12345678'),
        index=sceneNames)
