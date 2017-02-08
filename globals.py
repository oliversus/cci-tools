def init():
    global data_folder, figuresDir, latex_variables, NA1, NA2, SIB, main_folder, sceneTime, sceneTimes
    main_folder = "/cmsaf/esa_doku/ESA_Cloud_cci/publications/CC4CL_paper/"
    data_folder = main_folder + "data/"
    figuresDir = main_folder + "figures/"
    latex_variables = dict()
    sceneTimes = ["07221915", "07222058", "07270810", "01071228"]
    NA1 = sceneTimes[0]
    NA2 = sceneTimes[1]
    SIB = sceneTimes[2]
    TIM = sceneTimes[3]
    sceneTime = ""
