import torch
import cv2
import numpy as np
import pandas as pd
from torch import nn
import torch.nn.functional as F
from torch.autograd import Variable
from torchvision import transforms
import matplotlib.pyplot as plt

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

batch_size = 2**10 # 2**10

image_size = (100, 100, 3) # 100, 100, 3
encoding_size = 16 # 64

def to_img(x):
    return np.moveaxis(x.numpy() * 255, 0, -1).astype(np.uint8)

# Function to get the length-to-width ratio of the largest contour in an image
def getLWR(img):
    # Convert image to grayscale if it has more than one channel
    if img.shape[2] > 1:
        img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    # Add a constant border to the image
    img = cv2.copyMakeBorder(img, 10, 10, 10, 10, cv2.BORDER_CONSTANT, value=0)
    # Apply binary thresholding
    _, img = cv2.threshold(img, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    # Invert image if the background is white
    if img[0, 0] == 255:
        img = cv2.bitwise_not(img)
    # Find contours in the image
    contours, _ = cv2.findContours(img, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    # Get the largest contour by area
    contour = max(contours, key=cv2.contourArea)
    # Get dimensions of the minimum area rectangle enclosing the contour
    w, h = cv2.minAreaRect(contour)[1]
    # Return the length-to-width ratio
    return max(w, h) / min(w, h)

def getRedness(img):
  # image should be rgb
  data = img.reshape(-1, img.shape[-1])
  data = data[np.array([np.all(i != [0,0,0]) for i in data])]
  # redness is distance from (160,20,20)
  distances = np.linalg.norm(data - [160,20,20], axis = 1)
  distances = 255 - np.mean(distances)
  return distances

class VAE(nn.Module):
    def __init__(self):
        super(VAE, self).__init__()
        self.fc1 = nn.Linear(np.prod(image_size), 400)
        self.fc21 = nn.Linear(400, encoding_size)
        self.fc22 = nn.Linear(400, encoding_size)
        self.fc3 = nn.Linear(encoding_size, 400)
        self.fc4 = nn.Linear(400, np.prod(image_size))

    def encode(self, x):
        h1 = F.relu(self.fc1(x))
        return self.fc21(h1), self.fc22(h1)

    def reparametrize(self, mu, logvar):
        std = logvar.mul(0.5).exp_()
        if torch.cuda.is_available():
            eps = torch.cuda.FloatTensor(std.size()).normal_()
        else:
            eps = torch.FloatTensor(std.size()).normal_()
        eps = Variable(eps)
        return eps.mul(std).add_(mu)

    def decode(self, z):
        h3 = F.relu(self.fc3(z))
        return F.sigmoid(self.fc4(h3))

    def forward(self, x):
        mu, logvar = self.encode(x)
        z = self.reparametrize(mu, logvar)
        return self.decode(z), mu, logvar

reconstruction_function = nn.MSELoss(reduction='sum')

def loss_function(recon_x, x, mu, logvar):
    """
    recon_x: generating images
    x: origin images
    mu: latent mean
    logvar: latent log variance
    """
    BCE = reconstruction_function(recon_x, x)  # loss
    # loss = 0.5 * sum(1 + log(sigma^2) - mu^2 - sigma^2)
    KLD_element = mu.pow(2).add_(logvar.exp()).mul_(-1).add_(1).add_(logvar)
    KLD = torch.sum(KLD_element).mul_(-0.5)
    # KL divergence
    return BCE + KLD

model = VAE()
model.load_state_dict(torch.load('1_Data/embedding_to_image_coefficients.pth', map_location=device))

def embedding_to_image(embedding):
    embedding = torch.tensor(embedding)
    embedding = embedding.to(next(model.parameters()).dtype)
    decoded_image = model.decode(embedding)
    decoded_image = decoded_image.view(1, -1)
    decoded_image = to_img(decoded_image.detach().cpu().view(image_size[::-1]))  # Move back to CPU for processing
    return decoded_image[:,:,::-1]

def show_image(img, title = ''):
    plt.imshow(img)
    plt.axis('off')
    if(title != ''):
        plt.title(title)
    plt.show()