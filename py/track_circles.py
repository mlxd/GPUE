import cv, numpy
img= cv.LoadImage("wfc_1000.png",cv.CV_LOAD_IMAGE_GRAYSCALE)
eig_image = cv.CreateImage(cv.GetSize(img), cv.IPL_DEPTH_32F, 1)
temp_image = cv.CreateImage(cv.GetSize(img), cv.IPL_DEPTH_32F, 1)

circles=cv.CreateMat(img.width,1,cv.CV_32FC3)
cv.HoughCircles(img,circles,cv.CV_HOUGH_GRADIENT,2,10, 200,100)
c=numpy.asarray(circles)
for (x) in c:
	print x
