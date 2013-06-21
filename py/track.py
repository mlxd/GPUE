import cv
img= cv.LoadImage("foo2.jpg",cv.CV_LOAD_IMAGE_GRAYSCALE)
eig_image = cv.CreateImage(cv.GetSize(img), cv.IPL_DEPTH_32F, 1)
temp_image = cv.CreateImage(cv.GetSize(img), cv.IPL_DEPTH_32F, 1)
for (x,y) in cv.GoodFeaturesToTrack(img, eig_image, temp_image, 300, 0.1, 1.0, None, 3, True):
  print "good feature at", x,y

