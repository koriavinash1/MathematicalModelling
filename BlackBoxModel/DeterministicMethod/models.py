import torch
import torch.nn as nn
import torch.nn.functional as F

class DnCNN(nn.Module):
    def __init__(self, channels, num_of_layers=17):
        super(DnCNN, self).__init__()
        kernel_size = 3
        padding = 1
        features = 64

        self.dncnn = nn.Sequential()
        layers.append(nn.Conv2d(in_channels=channels, out_channels=features, kernel_size=kernel_size, padding=padding, bias=False))
        layers.append(nn.ReLU(inplace=True))
        for _ in range(num_of_layers-2):
            layers.append(nn.Conv2d(in_channels=features, out_channels=features, kernel_size=kernel_size, padding=padding, bias=False))
            layers.append(nn.BatchNorm2d(features))
            layers.append(nn.ReLU(inplace=True))
        layers.append(nn.Conv2d(in_channels=features, out_channels=channels, kernel_size=kernel_size, padding=padding, bias=False))
    def forward(self, x):
        out = self.dncnn(x)
        return out


class _Layer(nn.Sequential):
    def __init__(self, num_input_features, num_output_features, pad = 0):
        super(_Layer, self).__init__()
        self.add_module('norm1', nn.BatchNorm2d(num_input_features)),
        self.add_module('relu1', nn.ReLU(inplace=True))
        self.add_module('conv1', nn.Conv2d(num_input_features, num_output_features, kernel_size=3, stride=1, bias=False, padding = pad)),

    def forward(self, x):
        new_features = super(_Layer, self).forward(x)
        return new_features


class DCNN(nn.Module):
    def __init__(self, in_channels= 3, num_of_layers = 17):
        super(DCNN, self).__init__()
        self.features     = nn.Sequential(OrderedDict([]))
        self.features.add_module('Layer_' + str(1), _Layer(in_channels, 64, pad = 1))

        for i in range(num_of_layers - 2):
            self.features.add_module('Layer_' + str(i + 2), _Layer(64, 64, pad = 1))

        self.features.add_module('Layer_' + str(num_of_layers), _Layer(64, in_channels, pad = 1))
        
    def forward(self, x):
        return self.features(x)