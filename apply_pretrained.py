# -*- coding: utf-8 -*-
import argparse
import torch.utils.data
import numpy as np
import random
import pickle
import matplotlib.pyplot as plt
from DscSTAR import training, save_decoder_output
from pre_processing import pre_processing

seed = 0
torch.manual_seed(seed)
prep = 'mc'
pre_process_paras = {'prep': prep}

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run scSTAR_K_model')
    parser.add_argument('--data-folder', default='./inputs/', help='Path to the data folder')
    parser.add_argument('--input-file1', default='case.csv', help='First input csv file name')
    parser.add_argument('--input-file2', default='ctr.csv', help='Second input csv file name')
    parser.add_argument('--output-folder', default='./outputs/applyed_model/', help='Path to the output folder')
    parser.add_argument('--output-name1', default='case.out.csv', help='First output csv file name')
    parser.add_argument('--output-name2', default='ctr.out.csv', help='Second output csv file name')
    parser.add_argument('--model-path', default='./outputs/trained_model.pth', help='Path to the trained model')
    args = parser.parse_args()

    if torch.cuda.is_available():
      model = torch.load(args.model_path)
    else:
      model = torch.load(args.model_path, map_location=torch.device('cpu'))

    plotpath = args.output_folder + "UMAPplot.png"
    data_folder = args.data_folder
    dataset_file_list = [args.input_file1, args.input_file2]
    dataset_file_list = [data_folder + f for f in dataset_file_list]

    dataset, muX = pre_processing(dataset_file_list, pre_process_paras)

    # output
    output_file1 = args.output_folder + args.output_name1
    output_file2 = args.output_folder + args.output_name2
    
    save_decoder_output(model, dataset, muX, output_file1, output_file2, plotpath)
    
