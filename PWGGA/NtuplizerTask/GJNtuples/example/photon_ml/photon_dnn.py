#!/usr/bin/env python

from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

# Deterministic behavior
import os
os.environ['TF_CUDNN_USE_AUTOTUNE'] = '0'
prng_seed = 1009732533  # OEIS A002205
import numpy
numpy.random.seed(prng_seed)
import tensorflow
tensorflow.set_random_seed(prng_seed)

batch_size = 144
nb_classes = 2
nb_epoch = 100

ndense = 50
dropout = 0.1
nlayer = 4

activation = 'relu'
alpha_leaky_relu = 0.1

nconv2d = 0
conv2d_size = (3, 3)

model_loss = 'categorical_crossentropy'
model_optimizer = 'rmsprop'

import sys, getopt

try:
    option, argument = getopt.getopt(
        sys.argv[1:], '',
        ['batch-size=', 'nb-epoch=', 'ndense=', 'dropout=',
         'nlayer=', 'loss=', 'optimizer=',
         'leaky-relu', 'prelu', 'elu', 'relu'])
except getopt.GetoptError:
    sys.exit()
print(option, file = sys.stderr)
for o, a in option:
    if o == '--batch-size':
        batch_size = int(a)
    elif o == '--nb-epoch':
        nb_epoch = int(a)
    elif o == '--ndense':
        ndense = int(a)
    elif o == '--dropout':
        dropout = float(a)
    elif o == '--nlayer':
        nlayer = int(a)
    elif o == '--loss':
        model_loss = a
    elif o == '--optimizer':
        model_optimizer = a
    elif o == '--leaky-relu':
        activation = 'leaky_relu'
    elif o == '--prelu':
        activation = 'prelu'
    elif o == '--elu':
        activation = 'elu'
    elif o == '--relu':
        activation = 'relu'

# Weights to depopulate nonprompt photons, such the E_T distribution
# matches that of the prompt ones. Fit is performed in ROOT using a
# spline

def reweight_prompt(inv_sqrt_x):
    param = (1.58972e+02, -2.21542e+02, 1.18650e+02, -2.72830e+01,
             2.24137e+00)
    x = inv_sqrt_x**(-2)
    r = math.exp(param[0] + param[1] * math.log(x) +
                 param[2] * math.log(x)**2 +
                 param[3] * math.log(x)**3 +
                 param[4] * math.log(x)**4)
    return 366.358683429 / r

def reweight_nonprompt(inv_sqrt_x):
    param = (-3.20883e-05, 1.68863e-03, -2.82423e-02, 2.39077e+01,
             8.93126e-01)
    x = inv_sqrt_x**(-2)
    r = 0.5 * (1 + math.erf(1e+4 * (param[3] - x))) * \
        (param[0] * x**4 + param[1] * x**3 + param[2] * x**2 -
         (4 * param[0] * (param[3])**3 +
          3 * param[1] * (param[3])**2 +
          2 * param[2] * (param[3])) * x +
         (param[4]) + 3 * param[0] * (param[3])**4 +
         2 * param[1] * (param[3])**3 + param[2] * (param[3])**2) + \
        0.5 * (1 + math.erf(1e+4 * (x - param[3]))) * param[4]
    r *= reweight_prompt(inv_sqrt_x)
    return r

def add_activation(activation):
    if activation.lower() == 'prelu':
        act = keras.layers.advanced_activations.PReLU(
            init = 'random_uniform')
        model.add(act)
    elif activation.lower() == 'leaky_relu':
        act = keras.layers.advanced_activations.LeakyReLU(
            alpha = alpha_leaky_relu)
        model.add(act)
    elif activation.lower() == 'elu':
        act = keras.layers.advanced_activations.ELU(
            alpha = alpha_elu)
        model.add(act)
    elif activation.lower() == 'relu':
        model.add(keras.layers.Activation('relu'))


import math, h5py

filename = '/mnt/by-uuid/ca988570-a1e6-4d46-a186-20f7c6704a5c/scratch/photon_ml_ppb_pythia_dpmjet_5x5.h5'

f = h5py.File(filename, 'r')

# A small range that contains many prompt and nonprompt photons
#X = f['X'][137000:137000+1000]
#y = f['y'][137000:137000+1000]

X = f['X'][:]
y = f['y'][:]

f.close()

# Apply weighting
column_inv_sqrt_pt = X.shape[1] - 4
print('column_inv_sqrt_pt =', column_inv_sqrt_pt, file = sys.stderr)
reweight_prompt = numpy.vectorize(reweight_prompt)
reweight_nonprompt = numpy.vectorize(reweight_nonprompt)

i = numpy.argwhere(
    numpy.any(
        numpy.concatenate(
            (numpy.expand_dims(numpy.all(
                numpy.concatenate(
                    ((y == 1.0).T,
                     numpy.expand_dims(
                         numpy.random.rand(X.shape[0]) <
                         reweight_prompt(X[:, column_inv_sqrt_pt]),
                         axis = 0))),
                axis = 0), axis = 0),
             numpy.expand_dims(numpy.all(
                 numpy.concatenate(
                     ((y == 0.0).T,
                      numpy.expand_dims(
                          numpy.random.rand(X.shape[0]) <
                          reweight_nonprompt(X[:, column_inv_sqrt_pt]),
                          axis = 0))),
                 axis = 0), axis = 0))), axis = 0)).flatten()

#X = X[i,:25]
#y = y[i]

# Split into training and test samples
X_train, y_train, X_test, y_test = X[0::2], y[0::2], X[1::2], y[1::2]

# Transform 1 vs. 2 photons into [0, 1]
y_train -= 1
y_test -= 1

nfeature = X_train.shape[1]

import keras.backend

if keras.backend.image_dim_ordering() == 'th':
    X_train = X_train.reshape(X_train.shape[0], nfeature)
    X_test = X_test.reshape(X_test.shape[0], nfeature)
else:
    X_train = X_train.reshape(X_train.shape[0], nfeature)
    X_test = X_test.reshape(X_test.shape[0], nfeature)

mean = X_train.mean(0)
X_train -= mean
X_test -= mean
std = X_train.std(0)
X_train /= std
X_test /= std

f = open('photon_discr_norm.h', 'w')
print('static const float mean[' + str(len(mean)) + '] = {' + ', '.join(map(lambda x: '%.8e' % x, mean)) + '};', file = f)
print('static const float std[' + str(len(std)) + '] = {' + ', '.join(map(lambda x: '%.8e' % x, std)) + '};', file = f)
f.close()
sys.exit(0)

X_train = X_train.astype('float32')
X_test = X_test.astype('float32')
print('X_train shape:', X_train.shape, file = sys.stderr)
print(X_train.shape[0], 'train samples', file = sys.stderr)
print(X_test.shape[0], 'test samples', file = sys.stderr)

import keras.models
import keras.layers
if activation != 'relu':
    import keras.layers.advanced_activations
import keras.utils

# convert class vectors to binary class matrices
Y_train = keras.utils.np_utils.to_categorical(y_train, nb_classes)
Y_test = keras.utils.np_utils.to_categorical(y_test, nb_classes)

model = keras.models.Sequential()

if nconv2d > 0:
    model.add(keras.layers.Reshape((5, 5, 1), input_shape = (nfeature,)))
    for i in range(nconv2d):
        model.add(keras.layers.Conv2D(ndense, conv2d_size))
        add_activation(activation)
        model.add(keras.layers.Dropout(dropout))
    model.add(keras.layers.Flatten())
else:
    model.add(keras.layers.Dense(ndense, input_shape = (nfeature,)))
model.add(keras.layers.Dense(ndense))
add_activation(activation)
model.add(keras.layers.Dropout(dropout))
for i in range(nlayer):
    model.add(keras.layers.Dense(ndense))
    add_activation(activation)
    model.add(keras.layers.Dropout(dropout))
model.add(keras.layers.Dense(2))
model.add(keras.layers.Activation('softmax'))

model.compile(loss = model_loss,
              optimizer = model_optimizer,
              metrics = ['accuracy'])

model.fit(X_train, Y_train,
          batch_size = batch_size, nb_epoch = nb_epoch, verbose = 1,
          validation_data = (X_test, Y_test))
score = model.evaluate(X_test, Y_test, verbose=0)
print('Test score:', score[0], file = sys.stderr)
print('Test accuracy:', score[1], file = sys.stderr)

model.save('photon_discr.h5')

import kerasify

kerasify.export_model(model, 'photon_discr.model')
