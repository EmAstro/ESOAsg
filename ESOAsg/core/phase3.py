import os
import requests

from ESOAsg import msgs
from ESOAsg.ancillary import checks
from ESOAsg import default


def get_header_from_archive(file_id, text_file=None, verbose=False):  # written by Ema. 04.03.2020
    r"""Download header for the ESO given a `file_id`

    Args:
        file_id (str): ESO file ID for which the header will be downloaded
        text_file (str): text file where the header will be downloaded. If `None` it will it will be set to t
            he same string `file_id` but with a `.hdr` extension
        verbose (bool): if set to `True` additional info will be displayed

    Returns:
        None

    """

    # checks for connection to ESO archive
    archive_url = default.get_value('eso_archive_url')
    if not checks.connection_to_website(archive_url, timeout=10):
        msgs.error('Cannot connect to the ESO archive website:\n {}'.format(archive_url))

    # checks for file id
    assert isinstance(file_id, list) or isinstance(file_id, str), 'file_id needs to be a str or a list'
    if isinstance(file_id, str):
        list_of_files = [file_id]
    else:
        list_of_files = file_id
    list_of_files = [files if not files.endswith('.fits') else files.replace('.fits', '') for files in list_of_files]

    # checks for text_file
    assert isinstance(text_file, list) or isinstance(text_file, str) or \
           isinstance(text_file, (type(None), bytes)), 'text_file needs to be a str or a list'
    if isinstance(text_file, str):
        if len(list_of_files) == 1:
            list_of_outputs = [text_file]
        else:
            list_of_outputs = [output + text_file for output in list_of_files]
    elif isinstance(text_file, list):
        if len(list_of_files) == len(text_file):
            list_of_outputs = text_file
        else:
            list_of_outputs = [files + '.hdr' for files in list_of_files]
    else:
        list_of_outputs = [files + '.hdr' for files in list_of_files]

    # Downloading headers
    for file_name, file_out in zip(list_of_files, list_of_outputs):
        if os.path.isfile(file_out):
            if verbose:
                msgs.warning('Overwriting existing text file: {}'.format(file_out))
            os.remove(file_out)
        url_for_header = archive_url + 'hdr?DpId=' + file_name
        response_url = requests.get(url_for_header, allow_redirects=True)
        # Removing html from text
        header_txt = response_url.text.split('<pre>')[1].split('</pre>')[0]
        if not header_txt.startswith('No info found for'):
            file_header = open(file_out, 'w')
            for line in header_txt.splitlines():
                file_header.write(line + '\n')
            file_header.close()
            if verbose:
                msgs.info('Header successfully saved in: {}'.format(file_out))
        else:
            msgs.warning('{} is not present in the ESO archive'.format(file_name))
    return
