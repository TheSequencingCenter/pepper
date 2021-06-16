from os.path import isfile, join
from os import listdir
from torch.utils.data import Dataset
from pepper_variant.modules.python.Options import ImageSizeOptions
import torch
import torchvision.transforms as transforms
import h5py
import gzip
import pickle
import sys
import numpy as np


def get_file_paths_from_directory(directory_path):
    """
    Returns all paths of files given a directory path
    :param directory_path: Path to the directory
    :return: A list of paths of files
    """
    file_paths = [join(directory_path, file) for file in listdir(directory_path) if isfile(join(directory_path, file))
                  and file[-6:] == 'pkl.gz']
    return file_paths


class SequenceDataset(Dataset):
    """
    Arguments:
        A HDF5 file path
    """
    def __init__(self, image_directory):
        self.transform = transforms.Compose([transforms.ToTensor()])

        all_records = []
        hdf_files = get_file_paths_from_directory(image_directory)

        self.total_records = 0
        for hdf5_file_path in hdf_files:
            with h5py.File(hdf5_file_path, 'r') as hdf5_file:
                if 'summaries' in hdf5_file:
                    region_names = list(hdf5_file['summaries'].keys())

                    for region_name in region_names:
                        image_shape = hdf5_file['summaries'][region_name]['images'].shape[0]

                        for image_index in range(0, image_shape):
                            all_records.append((hdf5_file_path, region_name, image_index))

                        self.total_records += image_shape

        self.all_images = all_records

    def __getitem__(self, index):
        # load the image
        hdf5_filepath, region_name, indx = self.all_images[index]

        with h5py.File(hdf5_filepath, 'r') as hdf5_file:
            image = hdf5_file['summaries'][region_name]['images'][indx][()]
            type_label = hdf5_file['summaries'][region_name]['type_labels'][indx][()]
            base_label = hdf5_file['summaries'][region_name]['base_labels'][indx][()]

        return image, base_label, type_label

    def __len__(self):
        return len(self.all_images)


class SequenceDatasetFake(Dataset):
    """
    Arguments:
        A HDF5 file path
    """
    def __init__(self, image_directory):
        # self.transform = transforms.Compose([transforms.ToTensor()])
        # self.transform = transforms.Compose([])
        pickle_files = get_file_paths_from_directory(image_directory)
        self.file_lengths = []
        total_images = 0

        for pickle_file in pickle_files:

            file_sub_image_lengths = []
            with gzip.open(pickle_file, "rb") as image_file:
                while True:
                    try:
                        images = pickle.load(image_file)
                        image_length = len(images)
                        # as it is a 0-based
                        file_sub_image_lengths.append(image_length - 1)
                        total_images += image_length
                    except EOFError:
                        break

                self.file_lengths.append(file_sub_image_lengths)

        self.total_images = total_images
        self.pickle_files = pickle_files

    @staticmethod
    def my_collate(batch):
        candidate = [item[0] for item in batch]
        images = [item[1] for item in batch]
        base_predictions = [item[2] for item in batch]
        type_predictions = [item[3] for item in batch]
        images = torch.LongTensor(images)

        return [candidate, images, base_predictions, type_predictions]


    def __getitem__(self, index):
        # load the image
        file_index = 0
        sub_file_index = 0
        image_index = index

        found_index = False
        for i, sub_file_lengths in enumerate(self.file_lengths):
            file_index = i
            for j, sub_file_length in enumerate(sub_file_lengths):
                sub_file_index = j
                if image_index > sub_file_lengths[sub_file_index]:
                    image_index -= sub_file_lengths[sub_file_index]
                else:
                    found_index = True
                    break

            if found_index:
                break

        image_file = gzip.open(self.pickle_files[file_index], "rb")
        candidates = []

        for i in range(0, sub_file_index + 1):
            candidates = pickle.load(image_file)

        candidate = candidates[image_index]
        image_file.close()
        contig = candidate.contig
        image = np.array(candidate.image_matrix)

        position = candidate.position
        type_label = candidate.type_label
        base_label = candidate.base_label
        candidates = candidate.candidates
        candidate_frequency = candidate.candidate_frequency
        depth = candidate.depth

        base_predictions = np.zeros(ImageSizeOptions.TOTAL_LABELS)
        base_predictions[base_label] = 1

        type_predictions = np.zeros(ImageSizeOptions.TOTAL_TYPE_LABELS)
        type_predictions[type_label] = 1
        # print("######")
        # print(contig)
        #
        # print(position)
        # print(depth)
        # print(candidates)
        # print(candidate_frequency)
        # print(type(image))
        # print(base_predictions)
        # print(type_predictions)
        # print("######")
        # exit(0)

        return candidate, image, base_predictions, type_predictions
        # return contig, depth, candidates, candidate_frequency, image, position, base_predictions, type_predictions

    def __len__(self):
        return self.total_images


class SequenceDatasetHP(Dataset):
    """
    Arguments:
        A HDF5 file path
    """
    def __init__(self, image_directory):
        self.transform = transforms.Compose([transforms.ToTensor()])
        file_image_pair = []

        hdf_files = get_file_paths_from_directory(image_directory)

        for hdf5_file_path in hdf_files:
            with h5py.File(hdf5_file_path, 'r') as hdf5_file:
                if 'summaries' in hdf5_file:
                    image_names = list(hdf5_file['summaries'].keys())

                    for image_name in image_names:
                        # for hp_tag 1
                        file_image_pair.append((hdf5_file_path, image_name, 1))
                        # for hp_tag 2
                        file_image_pair.append((hdf5_file_path, image_name, 2))
                else:
                    sys.stderr.write("WARN: NO IMAGES FOUND IN FILE: " + hdf5_file_path + "\n")

        self.all_images = file_image_pair

    def __getitem__(self, index):
        # load the image
        hdf5_filepath, image_name, hp_tag = self.all_images[index]

        with h5py.File(hdf5_filepath, 'r') as hdf5_file:
            if hp_tag == 1:
                image = hdf5_file['summaries'][image_name]['image_hp1'][()]
                label = hdf5_file['summaries'][image_name]['label_hp1'][()]
            else:
                image = hdf5_file['summaries'][image_name]['image_hp2'][()]
                label = hdf5_file['summaries'][image_name]['label_hp2'][()]

        return image, label

    def __len__(self):
        return len(self.all_images)