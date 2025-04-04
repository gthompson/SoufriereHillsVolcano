# -*-coding:Utf-8 -*

# Copyright: Marielle MALFANTE - GIPSA-Lab -
# Univ. Grenoble Alpes, CNRS, Grenoble INP, GIPSA-lab, 38000 Grenoble, France
# (04/2018)
#
# marielle.malfante@gipsa-lab.fr (@gmail.com)
#
# This software is a computer program whose purpose is to automatically
# processing time series (automatic classification, detection). The architecture
# is based on machine learning tools.
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL

import json
from os.path import isfile, isdir
import datetime
from features import FeatureVector
import pickle
import numpy as np
from DataReadingFunctions import requestObservation, read_montserrat
from sklearn import preprocessing
from tools import butter_bandpass_filter
from featuresFunctions import energy, energy_u
from math import sqrt
from sklearn.metrics import confusion_matrix, accuracy_score
from sklearn.model_selection import cross_val_score, StratifiedShuffleSplit
from tools import print_cm
import time
from tools import extract_features
from copy import deepcopy

import os
import pandas as pd # added by Glenn
import glob


class Analyzer:
    """ Object containing the needed tools to analyze a Dataset.
    It contains the features scaler, the model, and the labels encoder,
    and can be used to train a model from supervised data.
    - scaler: None if learn has not been called, the learnt scaler otherwise
    - model: None if learn has not been called, the learnt model otherwise
    - labelEncoder: None if learn has not been called, the label encoder (label proper translation
    to int) otherwise
    - pathToCatalogue: path to the labeling catalogue. /!\ Catalogue should have a definite format.
    Check out README for more information on the catalogue shape.
    - catalogue: loaded catalogue of labels
    - _verbatim: how chatty do you want your Analyze to be?
    """

    def __init__(self, config, verbatim=0, catalog=None): #catalog added by Glenn
        """
        Initialization method
        """
        self.scaler = None
        self.model = deepcopy(config.learning['algo'])
        self.labelEncoder = None
        self.pathToCatalogue = config.general['project_root']+config.application['name'].upper()+'/'+config.learning['path_to_catalogue']
        print(self.pathToCatalogue)
        
        if isinstance(catalog, pd.DataFrame): # added by Glenn
            self.catalogue = catalog
        else:
            self.catalogue = pickle.load(open(self.pathToCatalogue,'rb'))
            #self.catalogue = pd.read_csv(open(self.pathToCatalogue,'rb'))
        self._verbatim = verbatim
        if self._verbatim>0:
            print('\n\n *** ANALYZER ***')
        return

    def __repr__(self):
        """
        Representation method (transform the object to str for display)
        """
        s = 'Analyzer object with model and scaler being: '+str(self.model)+' and ' +str(self.scaler)
        s += '\nCatalogue is at %s'%self.pathToCatalogue
        return s

    def learn(self, config, verbatim=None, forModelSelection=False, sss=None, model=None, featuresIndexes=None, returnData=False):
        """
        Method to train the analyzer.
        Labeled data are read from the catalogue, the data are preprocessed,
        features are extracted and scaled, and model is trained (with the
        stardard labels).
        All the arguments with default values are for a "non classic" use of
        the analyzer object (model selection for example)
        Return None, but can return the data and labels if specified in returnData.
        """
        if verbatim is None:
            verbatim=self._verbatim

        # Get or define usefull stuff
        features = FeatureVector(config, verbatim=verbatim)
        nData = len(self.catalogue.index)
        if returnData:
            allData = np.zeros((nData,),dtype=object)
        allLabels = np.zeros((nData,),dtype=object)
        allFeatures = np.zeros((nData,features.n_domains*features.n_features),dtype=float)

        # Read all labeled signatures (labels+data) from the catalogue, and extract features
        tStart = time.time()
        catalog_length = len(self.catalogue.index)
        WAVTOPDIR = config.data_to_analyze['path_to_data'] # Glenn. path has miniseed_c hardcoded at start. I want to change this to whatever the config says
        for i in range(catalog_length):
            if verbatim > 1:
                print('Processing waveform %d of %d' % (i, catalog_length))
            if verbatim > 3:
                print(self.catalogue.iloc[i].to_string())
            secondFloat = self.catalogue.iloc[i]['second']
            tStartSignature = datetime.datetime(int(self.catalogue.iloc[i]['year']),     \
                                                int(self.catalogue.iloc[i]['month']),    \
                                                int(self.catalogue.iloc[i]['day']),      \
                                                int(self.catalogue.iloc[i]['hour']),     \
                                                int(self.catalogue.iloc[i]['minute']),   \
                                                int(secondFloat), \
                                                int((secondFloat-int(secondFloat))*1000000)) #microseconds
            duration = self.catalogue.iloc[i]['length']
            path = self.catalogue.iloc[i]['path']     
            path = path.replace('miniseed_c', WAVTOPDIR)
            
            # Get label and check that it is single label (multi label not supported yet)
            lab = self.catalogue.iloc[i]['class']
            if type(lab) is list:
                print('Multi label not implemented for learning yet')
                return None
            allLabels[i] = lab
            
            # LOAD WAVEFORM
            #(fs, signature) = requestObservation(config, tStartSignature, duration, path, verbatim=0)
            # first check for features file
            
            # what traceID are we looking for - should figure out how to write this into the config
            fptr = open('./AAA-master/MONTSERRAT/current_traceID.txt','r')
            traceID = fptr.read()
            fptr.close()     
            
            if traceID=='*':
                mseedbase = os.path.basename(path)
                allpkls = glob.glob(os.path.join('features',mseedbase.replace('.mseed', '*.%s.pkl')))
                for featurespkl in allpkls:
                    print('Reading %s' % featurespkl)
                    featuresList = pd.read_pickle(featurespkl)
            
                    # APPEND FEATURES FOR THIS EVENT+TRACE TO ALLFEATURES
                    allFeatures[i] = featuresList[0] # subscript 0 not needed
                    if verbatim > 2:
                        print('Got %d features for %s in %s' % (len(featuresList), traceID, path))  
                    if verbatim > 3:
                        print(featuresList)
        
        # OUTPUT COMPUTATION TIME
        tEnd = time.time()
        if verbatim>0:
            print('Training data have been read and features have been extracted ', np.shape(allFeatures))
            print('Computation time for reading and feature computation: ', tEnd-tStart)


        # Compress labels and features in case of None values (if reading is empty for example)
        i = np.where(allLabels != np.array(None))[0]
        allFeatures = allFeatures[i]
        allLabels = allLabels[i]
        if returnData:
            allData = allData[i]

        # Transform labels
        self.labelEncoder = preprocessing.LabelEncoder().fit(allLabels)
        allLabelsStd = self.labelEncoder.transform(allLabels)
        if verbatim>0:
            print('Model will be trained on %d classes'%len(self.labelEncoder.classes_), np.unique(allLabelsStd), self.labelEncoder.classes_)

        # Scale features and store scaler
        self.scaler = preprocessing.StandardScaler().fit(allFeatures)
        allFeatures = self.scaler.transform(allFeatures)
        if verbatim>0:
            print('Features have been scaled')

        # Get model from learning configuration file and learn
        self.model = deepcopy(config.learning['algo'])

        if forModelSelection:
            if model is None:
                pass
            else:
                self.model = model

        tStartLearning = time.time()
        if featuresIndexes is None:
            self.model = self.model.fit(allFeatures, allLabelsStd)
        else:
            self.model = self.model.fit(allFeatures[:,featuresIndexes], allLabelsStd)
        tEndLearning = time.time()

        #  Model Evaluation (a) with score, (b) with X-validation
        # NB: When model is trained (and evaluated by X-validation or score),
        # threshold is NOT used. Threshold is only used when the 'unknown'
        # class can occur (and this is obvisouly not the case with supervised
        # training)
        print('Model has been trained: ', self.model)
        print('Computation time: ', tEndLearning-tStartLearning)

        if featuresIndexes is None:
            allPredictions = self.model.predict(allFeatures)
            allProbabilities = self.model.predict_proba(allFeatures) # added by Glenn
        else:
            allPredictions = self.model.predict(allFeatures[:,featuresIndexes])
            allProbabilities = self.model.predict_proba(allFeatures[:,featuresIndexes])  # added by Glenn

        # (a) Score evaluation
        print('Model score is: ', accuracy_score(allLabelsStd,allPredictions))
        lab = list(range(len(self.labelEncoder.classes_))) # 'unknown' class not needed.
        CM = confusion_matrix(allLabelsStd,allPredictions,labels=lab)
        print('and associated confusion matrix is:')
        print_cm(CM, list(self.labelEncoder.classes_),hide_zeroes=True,max_str_label_size=2,float_display=False)

        # (b) X-validation
        sss = config.learning['cv']
        print(sss)
        CM=list()
        acc=list()
        model_Xval = deepcopy(self.model)
        for (i, (train_index, test_index)) in enumerate(sss.split(allFeatures, allLabelsStd)):
            predictionsStd = model_Xval.fit(allFeatures[train_index], allLabelsStd[train_index]).predict(allFeatures[test_index])
            predictions = self.labelEncoder.inverse_transform(predictionsStd)
            CM.append(confusion_matrix(allLabels[test_index],predictions, labels=self.labelEncoder.classes_))
            acc.append(accuracy_score(allLabels[test_index],predictions))
        print('Cross-validation results: ', np.mean(acc)*100, ' +/- ', np.std(acc)*100, ' %')
        print_cm(np.mean(CM, axis=0),self.labelEncoder.classes_,hide_zeroes=True,max_str_label_size=2,float_display=False)

        if returnData:
            return allData, allLabels, acc, allPredictions, allProbabilities
        else:
            return None

    def save(self, config):
        """
        Method used to save the object for later use (depending on the
        application, training can take a while and you might want to save the analyzer)
        """
        path = config.general['project_root'] + config.application['name'].upper() + '/res/' + config.configuration_number + '/' + config.general['path_to_res']
        savingPath = path+'analyzer'
        pickle.dump(self.__dict__,open(savingPath,'wb'),2)
        if self._verbatim > 0:
            print('Analyzer has been saved at: ', savingPath)
        return

    def load(self, config):
        """
        Method used to load the object.
        """
        verbatim = self._verbatim
        path = config.general['project_root'] + config.application['name'].upper() + '/res/' + config.configuration_number + '/' + config.general['path_to_res']
        savingPath = path+'analyzer'
        tmp_dict = pickle.load(open(savingPath,'rb'))
        self.__dict__.update(tmp_dict)
        self._verbatim = verbatim
        if self._verbatim > 0:
            print('Analyzer has been loaded from: ', savingPath)
        return
