#!/bin/bash

delLatLon=0.1

# delLat delLon doRGB sceneTime corrected plotCot plotCalipso

./paperFigures.py ${delLatLon} ${delLatLon} True 07222058 False False True
./paperFigures.py ${delLatLon} ${delLatLon} True 07221915 False False True
./paperFigures.py ${delLatLon} ${delLatLon} True 07270810 False False True

# # case 07222058
# # paper plot 1: corrected Ctp, no Cot
# ./paperFigures.py ${delLatLon} ${delLatLon} False 07222058 True False True
# # paper plot 2: uncorrected Ctp, no Cot
# ./paperFigures.py ${delLatLon} ${delLatLon} False 07222058 False False True
# # internal plot 1: corrected Ctp, Cot
# ./paperFigures.py ${delLatLon} ${delLatLon} False 07222058 True True True
# # internal plot 2: uncorrected Ctp, Cot
# ./paperFigures.py ${delLatLon} ${delLatLon} False 07222058 False True True

# # case 07270810
# # paper plot 1: corrected Ctp, no Cot
# ./paperFigures.py ${delLatLon} ${delLatLon} False 07270810 True False True
# # paper plot 2: uncorrected Ctp, no Cot
# ./paperFigures.py ${delLatLon} ${delLatLon} False 07270810 False False True
# # internal plot 1: corrected Ctp, Cot
# ./paperFigures.py ${delLatLon} ${delLatLon} False 07270810 True True True
# # internal plot 2: uncorrected Ctp, Cot
# ./paperFigures.py ${delLatLon} ${delLatLon} False 07270810 False True True

# # case 07221915
# # paper plot 1: corrected Ctp, no Cot
# ./paperFigures.py ${delLatLon} ${delLatLon} False 07221915 True False True
# # paper plot 2: uncorrected Ctp, no Cot
# ./paperFigures.py ${delLatLon} ${delLatLon} False 07221915 False False True
# # internal plot 1: corrected Ctp, Cot
# ./paperFigures.py ${delLatLon} ${delLatLon} False 07221915 True True True
# # internal plot 2: uncorrected Ctp, Cot
# ./paperFigures.py ${delLatLon} ${delLatLon} False 07221915 False True True


