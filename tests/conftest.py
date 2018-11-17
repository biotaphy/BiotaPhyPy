"""Test configuration fixtures
"""
import glob
import os


import pytest

# .............................................................................
# .                                 Constants                                 .
# .............................................................................
ALIGNMENTS_DIR = 'alignments'
ANC_STATE_PACKAGES_DIR = 'ancestral_state_packages'
ANC_DIST_PACKAGES_DIR = 'ancestral_distribution_packages'
TREES_DIR = 'trees'

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
SAMPLE_DATA_PATH = os.path.join(THIS_DIR, 'data_dir')


# .............................................................................
class SampleDataFiles(object):
    """This class is used to retrieve sample data for the tests.

    Note:
        * For test files, the format should be something like:
            "(in)valid_{name}.{extension}"
    """
    # .....................................
    def get_alignments(self, fmt, is_valid):
        """Gets an alignment file from the sample data.

        Args:
            fmt :  The format of the file you want (csv, json, phylip, table)
            is_valid : Return valid data files if true, invalid files if not

        Returns:
            A list of alignments matching the arguments
        """
        ALIGNMENTS_PATH = os.path.join(SAMPLE_DATA_PATH, ALIGNMENTS_DIR)
        return glob.iglob(
            self._get_glob_string(ALIGNMENTS_PATH, is_valid,
                                  self._get_format_extension(fmt)))

    # .....................................
    def get_ancestral_distribution_packages(self, is_valid):
        """Gets a list of the available packages for testing

        Args:
            is_valid : Return valid data files if true, invalid files if not

        Returns:
            A list of (tree, alignment, results) tuples
        """
        packages_path = os.path.join(SAMPLE_DATA_PATH, ANC_DIST_PACKAGES_DIR)
        return self._get_packages(packages_path, is_valid)

    # .....................................
    def get_ancestral_state_packages(self, is_valid):
        """Gets a list of the available packages for testing

        Args:
            is_valid : Return valid data files if true, invalid files if not

        Returns:
            A list of (tree, alignment, results) tuples
        """
        packages_path = os.path.join(SAMPLE_DATA_PATH, ANC_STATE_PACKAGES_DIR)
        return self._get_packages(packages_path, is_valid)

    # .....................................
    def get_trees(self, fmt, is_valid):
        """Gets an alignment file from the sample data.

        Args:
            fmt : The format of the file you want (newick, nexus)
            is_valid : Return valid data files if true, invalid files if not

        Returns:
            A list of tree filenames
        """
        TREE_PATH = os.path.join(SAMPLE_DATA_PATH, TREES_DIR)
        return glob.iglob(
            self._get_glob_string(TREE_PATH, is_valid,
                                  self._get_format_extension(fmt)))

    # .....................................
    def _get_format_extension(self, fmt):
        if fmt.lower() == 'csv':
            return '.csv'
        elif fmt.lower() == 'json':
            return '.json'
        elif fmt.lower() == 'newick':
            return '.tre'
        elif fmt.lower() == 'nexus':
            return '.nex'
        elif fmt.lower() == 'phylip':
            return '.phylip'
        elif fmt.lower() == 'table':
            return '.tbl'
        else:
            raise Exception('Cannot handle format: {}'.format(fmt))

    # .....................................
    def _get_glob_string(self, search_dir, is_valid, fmt_ext):
        if is_valid:
            valid_str = 'valid'
        else:
            valid_str = 'invalid'
        return os.path.join(search_dir, '{}_*{}'.format(valid_str, fmt_ext))

    # .....................................
    def _get_packages(self, packages_path, is_valid):
        """Gets a list of the available packages for testing

        Args:
            packages_path : A directory to look for packages within
            is_valid : Return valid data files if true, invalid files if not


        Returns:
            A list of (tree, alignment, results) tuples
        """
        packages = []
        package_dirs = glob.glob(
            self._get_glob_string(packages_path, is_valid, ''))
        for pkg_dir in package_dirs:
            tree_fn = None
            align_fn = None
            results_fn = None
            for fn in glob.glob(os.path.join(pkg_dir, '*')):
                basename = os.path.basename(fn)
                if basename.lower().startswith('tree'):
                    tree_fn = fn
                elif basename.lower().startswith('align'):
                    align_fn = fn
                elif basename.lower().startswith('result'):
                    results_fn = fn
            if tree_fn is not None and align_fn is not None and \
                    results_fn is not None:
                packages.append((tree_fn, align_fn, results_fn))
        return packages


# .............................................................................
@pytest.fixture(scope="session")
def data_files():
    """Gets test fixture used to retrieve sample data files.

    Returns:
        A `SampleDataFiles` object
    """
    return SampleDataFiles()
