import helpers as hlp
from sys import platform
import os


def make_read_lists(folder_path, load_filepaths_from_text_file=False):
    """
    This function collects file paths for the Teloprime cDNA and
    256 Cell and Shilesd RNA data set.

    :param load_filepaths_from_text_file: If True then loads files from a pre-saved text file
    """

    if platform == "linux" or platform == "linux2":
        mac_or_linux = 'linux'
    elif platform == "darwin":
        mac_or_linux = 'mac'
    else:
        mac_or_linux = 'windows'
        print("Check file paths if they are proper")

    if mac_or_linux == 'mac':
        data_path = os.path.join("/Users/adnaniazi/mnt/kjempetuja", folder_path) 

    else:
        data_path = folder_path 


    #-------------------------------------------------------------------------------------------------------------------
    cwd_directory_path = os.path.join(os.getcwd(), 'file_paths')

    if load_filepaths_from_text_file:
        # Return file paths saved in a text file as a list
        data_path = hlp.load_filepaths_from_json('files_paths.json', directory_path=cwd_directory_path)

    else:
        files_paths = hlp.find_all_fast5s(data_path)

        # Save these paths to JSON files in a file_paths directory
        if not os.path.exists(cwd_directory_path):
            os.makedirs(cwd_directory_path)

        hlp.save_filepaths_in_json(files_paths,
                                   directory_path=cwd_directory_path,
                                   file_name='files_paths.json')

    return files_paths



