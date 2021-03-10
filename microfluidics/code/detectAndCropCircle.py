import cv2
import numpy as np
import sys
import argparse
import os
import glob

# construct the argument parser and parse the arguments
ap = argparse.ArgumentParser()
#ap.add_argument("-i", "--image", required = True, help = "Path to the image")
ap.add_argument("-i", "--input", required = True, help = "Directory to images")
ap.add_argument("-o", "--output",required = True, help = "Directory to write cropped images")
args = vars(ap.parse_args())


#Check if output directory exists to prevent overwriting
try:
    os.makedirs(args["output"])
except OSError as e:
    if os.path.exists(args["output"]):
        sys.exit("Error creating output directory; already exists")


#img1 = cv2.imread(args["image"])
#img = cv2.imread(args["image"],0)
#print(os.listdir(args["input"]))


images=[]
for filename in os.listdir(args["input"]):
    #print(args["input"]+"/"+filename)
    images.append(args["input"]+"/"+filename)

#print(images)

for image in images:
    print(image)
    img1 = cv2.imread(image)
    img = cv2.imread(image,0)

    gray = cv2.cvtColor(img1, cv2.COLOR_BGR2GRAY)
    ret, thresh = cv2.threshold(gray, 50, 255, cv2.THRESH_BINARY)

    # Create mask
    height,width = img.shape
    mask = np.zeros((height,width), np.uint8)

    edges = cv2.Canny(thresh, 100, 200)
    #cv2.imshow('detected ',gray)
    cimg=cv2.cvtColor(img, cv2.COLOR_GRAY2BGR)
    #circles = cv2.HoughCircles(edges, cv2.HOUGH_GRADIENT, 1, 10000, param1 = 50, param2 = 30, minRadius = 0, maxRadius = 0)
    circles = cv2.HoughCircles(gray, cv2.HOUGH_GRADIENT, 1.2, 100)
    if(circles is not None):
        circles = np.round(circles[0, :]).astype("int")
        ##CHANGE if you want all circles in image or just one
        #for (x, y, r) in circles:
        #    cv2.circle(mask,(x,y),r,(255,255,255),thickness=-1)
        #This option for just one circle from image
        (x,y,r) = circles[0]
        cv2.circle(mask,(x,y),r,(255,255,255),thickness=-1)

        # Copy that image using that mask
        masked_data = cv2.bitwise_and(img1, img1, mask=mask)

        # Apply Threshold
        _,thresh = cv2.threshold(mask,1,255,cv2.THRESH_BINARY)

        # Find Contour
        contours, _ = cv2.findContours(thresh,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_SIMPLE)

        x,y,w,h = cv2.boundingRect(contours[0])

        # Crop masked_data
        crop = masked_data[y:y+h,x:x+w]

        #Code to close Window
        #cv2.imshow('detected Edge',img1)
        #cv2.imshow('Cropped Eye',crop)

        #cv2.waitKey(0)
        #cv2.destroyAllWindows()
        #filename = args["image"]
        filepath = os.path.splitext(image)[0]
        #filepath = os.path.splitext(args["image"])[0]
        filename = os.path.basename(filepath)
        cv2.imwrite(args["output"]+'/'+filename+"_crop.tif",crop)
