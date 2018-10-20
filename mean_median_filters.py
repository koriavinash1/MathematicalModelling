
# coding: utf-8

# In[22]:


import cv2
import numpy as np
import matplotlib.pyplot as plt

noise_level = 150.
image = cv2.imread('./images/08.png', 0)
noise = np.random.normal(0, noise_level/5., image.shape)
image_noise = image + noise

kernel = np.ones((5,5),np.float32)/25.0
processed_image_mean = cv2.filter2D(np.uint8(image_noise),-1,kernel)

processed_image_median = cv2.medianBlur(np.uint8(image_noise), 3)

plt.figure(figsize = (15, 15))
plt.subplot(1, 3, 1)
plt.title("Noisy Image")
plt.imshow(image_noise)

plt.subplot(1, 3, 2)
plt.title("Median Filter Processed Image")
plt.imshow(processed_image_median)

plt.subplot(1, 3, 3)
plt.title("Mean Filter Processed Image")
plt.imshow(processed_image_mean)

plt.show()

