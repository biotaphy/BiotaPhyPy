"""This module contains classes and functions for testing data readers

Notes:
    * Uses pytest style testing
"""
try:
    from StringIO import StringIO
    string_formats = (basestring)
except:
    # Python 3
    from io import StringIO
    string_formats = (str)

import pytest

import ancestral_reconstruction.helpers.data_readers as dr
from ancestral_reconstruction.helpers.sequence import Sequence


# .............................................................................
class Test_create_sequence_list_from_dict(object):
    """Test class for the create_sequence_list_from_dict method
    """
    # .....................................
    def test_invalid_dict(self):
        """Tests that function fails with a dict that has the wrong structure
        """
        test_dict = {
            "name1": [1, 2, 3],
            "name2": "bad value"
        }
        with pytest.raises(dr.AlignmentIOError):
            dr.create_sequence_list_from_dict(test_dict)

    # .....................................
    def test_valid_dict(self):
        """Tests that the function operates correctly with valid input data
        """
        test_dict = {
            "A": [0.9, 0.2, 0.2, 0.3, 0.4, 0.4],
            "B": [0.01, 0.1, 0.2, 0.3, 0.4, 0.4],
            "C": [0.8, 0.1, 0.2, 0.3, 0.4, 0.4],
            "D": [0.3, 0.1, 0.2, 0.3, 0.4, 0.4],
            "E": [0.001, 0.1, 0.2, 0.3, 0.4, 0.4],
            "F": [0.11, 0.1, 0.2, 0.3, 0.4, 0.4],
            "G": [0.99, 0.2, 0.2, 0.3, 0.4, 0.4]
        }
        sequence_list, headers = dr.create_sequence_list_from_dict(test_dict)
        for i in sequence_list:
            assert isinstance(i, Sequence)
            assert isinstance(i.name, string_formats)
            assert len(i.cont_values) > 0


# .............................................................................
class Test_read_csv_alignment_flo(object):
    """Tests the dr.read_csv_alignment_flo method
    """
    # .....................................
    def test_file_invalid(self, data_files):
        """Tests that the invalid alignment files fail properly

        Args:
            data_files : A pytest fixture defined in conftest.py
        """
        invalid_csv_alignments = data_files.get_alignments('csv', False)
        for fn in invalid_csv_alignments:
            with open(fn) as in_csv:
                with pytest.raises(dr.AlignmentIOError):
                    dr.read_csv_alignment_flo(in_csv)

    # .....................................
    def test_file_valid(self, data_files):
        """Tests that the valid alignment files do not fail

        Args:
            data_files : A pytest fixture defined in conftest.py
        """
        valid_csv_alignments = data_files.get_alignments('csv', True)
        for fn in valid_csv_alignments:
            with open(fn) as in_csv:
                try:
                    sequence_list, headers = dr.read_csv_alignment_flo(in_csv)
                    assert len(headers) > 0
                    assert len(sequence_list) > 0
                    for i in sequence_list:
                        assert isinstance(i, Sequence)
                        assert isinstance(i.name, string_formats)
                        assert len(i.cont_values) > 0
                except Exception as e:
                    print('Raised exception: {}'.format(str(e)))
                    assert False

    # .....................................
    def test_stringio_invalid(self, data_files):
        """Tests that invalid alignment files cause the proper failure

        Tests that the invalid alignment files fail properly when loaded into
        StringIO objects

        Args:
            data_files : A pytest fixture defined in conftest.py
        """
        invalid_csv_alignments = data_files.get_alignments('csv', False)
        for fn in invalid_csv_alignments:
            with open(fn) as in_csv:
                csv_stringio = StringIO()
                csv_stringio.write(in_csv.read())
                csv_stringio.seek(0)
                with pytest.raises(dr.AlignmentIOError):
                    dr.read_csv_alignment_flo(csv_stringio)

    # .....................................
    def test_stringio_valid(self, data_files):
        """Tests that alignments are read properly with valid input data

        Tests that the valid alignment files do not fail when loaded into
        StringIO objects

        Args:
            data_files : A pytest fixture defined in conftest.py
        """
        valid_csv_alignments = data_files.get_alignments('csv', True)
        for fn in valid_csv_alignments:
            with open(fn) as in_csv:
                csv_stringio = StringIO()
                csv_stringio.write(in_csv.read())
                csv_stringio.seek(0)
                try:
                    (sequence_list, headers) = dr.read_csv_alignment_flo(
                        csv_stringio)
                    assert len(headers) > 0
                    assert len(sequence_list) > 0
                    for i in sequence_list:
                        assert isinstance(i, Sequence)
                        assert isinstance(i.name, string_formats)
                        assert len(i.cont_values) > 0
                except Exception as e:
                    print('Raised exception: {}'.format(str(e)))
                    assert False


# .............................................................................
class Test_read_json_alignment_flo(object):
    """Test class for the dr.read_json_alignment_flo method
    """
    # .....................................
    def test_file_invalid(self, data_files):
        """Tests that the invalid alignment files fail properly

        Args:
            data_files : A pytest fixture defined in conftest.py
        """
        invalid_json_alignments = data_files.get_alignments('json', False)
        for fn in invalid_json_alignments:
            with open(fn) as in_json:
                with pytest.raises(dr.AlignmentIOError):
                    dr.read_json_alignment_flo(in_json)

    # .....................................
    def test_file_valid(self, data_files):
        """Tests that the valid alignment files do not fail

        Args:
            data_files : A pytest fixture defined in conftest.py
        """
        valid_json_alignments = data_files.get_alignments('json', True)
        for fn in valid_json_alignments:
            with open(fn) as in_json:
                try:
                    (sequence_list, headers
                     ) = dr.read_json_alignment_flo(in_json)
                    if headers is not None:
                        assert len(headers) > 0
                    assert len(sequence_list) > 0
                    for i in sequence_list:
                        assert isinstance(i, Sequence)
                        assert isinstance(i.name, string_formats)
                        assert len(i.cont_values) > 0
                except Exception as e:
                    print('Raised exception: {}'.format(str(e)))
                    assert False

    # .....................................
    def test_stringio_invalid(self, data_files):
        """Tests that invalid alignment files fail properly

        Tests that the invalid alignment files fail properly when loaded into
        StringIO objects

        Args:
            data_files : A pytest fixture defined in conftest.py
        """
        invalid_json_alignments = data_files.get_alignments('json', False)
        for fn in invalid_json_alignments:
            with open(fn) as in_json:
                json_stringio = StringIO()
                json_stringio.write(in_json.read())
                json_stringio.seek(0)
                with pytest.raises(dr.AlignmentIOError):
                    dr.read_json_alignment_flo(json_stringio)

    # .....................................
    def test_stringio_valid(self, data_files):
        """Tests that valid JSON alignment files do not fail

        Tests that the valid alignment files do not fail when loaded into
        StringIO objects

        Args:
            data_files : A pytest fixture defined in conftest.py
        """
        valid_json_alignments = data_files.get_alignments('json', True)
        for fn in valid_json_alignments:
            with open(fn) as in_json:
                json_stringio = StringIO()
                json_stringio.write(in_json.read())
                json_stringio.seek(0)
                try:
                    (sequence_list, headers
                     ) = dr.read_json_alignment_flo(json_stringio)
                    if headers is not None:
                        assert len(headers) > 0
                    assert len(sequence_list) > 0
                    for i in sequence_list:
                        assert isinstance(i, Sequence)
                        assert isinstance(i.name, string_formats)
                        assert len(i.cont_values) > 0
                except Exception as e:
                    print('Raised exception: {}'.format(str(e)))
                    assert False


# .............................................................................
class Test_read_phylip_alignment_flo(object):
    """Test class for the dr.read_phylip_alignment_flo method
    """
    # .....................................
    def test_file_invalid(self, data_files):
        """Tests that the invalid alignment files fail properly

        Args:
            data_files : A pytest fixture defined in conftest.py
        """
        invalid_phylip_alignments = data_files.get_alignments('phylip', False)
        for fn in invalid_phylip_alignments:
            with open(fn) as in_phylip:
                with pytest.raises(dr.AlignmentIOError):
                    dr.read_phylip_alignment_flo(in_phylip)

    # .....................................
    def test_file_valid(self, data_files):
        """Tests that the valid alignment files do not fail

        Args:
            data_files : A pytest fixture defined in conftest.py
        """
        valid_phylip_alignments = data_files.get_alignments('phylip', True)
        for fn in valid_phylip_alignments:
            with open(fn) as in_phylip:
                try:
                    sequence_list = dr.read_phylip_alignment_flo(in_phylip)
                    assert len(sequence_list) > 0
                    for i in sequence_list:
                        assert isinstance(i, Sequence)
                        assert isinstance(i.name, string_formats)
                        assert i.seq is not None
                except Exception as e:
                    print('Raised exception: {}'.format(str(e)))
                    assert False

    # .....................................
    def test_stringio_invalid(self, data_files):
        """Test that invalid phylip alignment files cause a failure

        Tests that the invalid alignment files fail properly when loaded into
        StringIO objects

        Args:
            data_files : A pytest fixture defined in conftest.py
        """
        invalid_phylip_alignments = data_files.get_alignments('phylip', False)
        for fn in invalid_phylip_alignments:
            with open(fn) as in_phylip:
                phylip_stringio = StringIO()
                phylip_stringio.write(in_phylip.read())
                phylip_stringio.seek(0)
                with pytest.raises(dr.AlignmentIOError):
                    dr.read_phylip_alignment_flo(phylip_stringio)

    # .....................................
    def test_stringio_valid(self, data_files):
        """Tests that valid phylip alignment files are loaded properly

        Tests that the valid alignment files do not fail when loaded into
        StringIO objects

        Args:
            data_files : A pytest fixture defined in conftest.py
        """
        valid_phylip_alignments = data_files.get_alignments('phylip', True)
        for fn in valid_phylip_alignments:
            with open(fn) as in_phylip:
                phylip_stringio = StringIO()
                phylip_stringio.write(in_phylip.read())
                phylip_stringio.seek(0)
                try:
                    sequence_list = dr.read_phylip_alignment_flo(
                        phylip_stringio)
                    assert len(sequence_list) > 0
                    for i in sequence_list:
                        assert isinstance(i, Sequence)
                        assert isinstance(i.name, string_formats)
                        assert i.seq is not None
                except Exception as e:
                    print('Raised exception: {}'.format(str(e)))
                    assert False


# .............................................................................
class Test_read_phylip_continuous_alignment_flo(object):
    """Test class for the read_phylip_continuous_alignment_flo method
    """
    # .....................................
    def test_file_invalid(self, data_files):
        pass

    # .....................................
    def test_file_valid(self, data_files):
        pass

    # .....................................
    def test_stringio_invalid(self, data_files):
        pass

    # .....................................
    def test_stringio_valid(self, data_files):
        pass


# .............................................................................
class Test_read_table_continuous_alignment_flo(object):
    """Test class for the read_table_continuous_alignment_flo method
    """
    # .....................................
    def test_file_invalid(self, data_files):
        """Tests that the invalid alignment files fail properly

        Args:
            data_files : A pytest fixture defined in conftest.py
        """
        invalid_table_alignments = data_files.get_alignments('table', False)
        for fn in invalid_table_alignments:
            with open(fn) as in_table:
                with pytest.raises(dr.AlignmentIOError):
                    dr.read_table_alignment_flo(in_table)

    # .....................................
    def test_file_valid(self, data_files):
        """Tests that the valid alignment files do not fail

        Args:
            data_files : A pytest fixture defined in conftest.py
        """
        valid_table_alignments = data_files.get_alignments('table', True)
        for fn in valid_table_alignments:
            with open(fn) as in_table:
                try:
                    sequence_list = dr.read_table_alignment_flo(in_table)
                    assert len(sequence_list) > 0
                    for i in sequence_list:
                        assert isinstance(i, Sequence)
                        assert isinstance(i.name, string_formats)
                        assert len(i.cont_values) > 0
                except Exception as e:
                    print('Raised exception: {}'.format(str(e)))
                    assert False

    # .....................................
    def test_stringio_invalid(self, data_files):
        """Tests that invalid table alignment files cause a failure

        Tests that the invalid alignment files fail properly when loaded into
        StringIO objects

        Args:
            data_files : A pytest fixture defined in conftest.py
        """
        invalid_table_alignments = data_files.get_alignments('table', False)
        for fn in invalid_table_alignments:
            with open(fn) as in_table:
                table_stringio = StringIO()
                table_stringio.write(in_table.read())
                table_stringio.seek(0)
                with pytest.raises(dr.AlignmentIOError):
                    dr.read_table_alignment_flo(table_stringio)

    # .....................................
    def test_stringio_valid(self, data_files):
        """Tests that valid table alignment files are loaded correctly.

        Tests that the valid alignment files do not fail when loaded into
        StringIO objects

        Args:
            data_files : A pytest fixture defined in conftest.py
        """
        valid_table_alignments = data_files.get_alignments('table', True)
        for fn in valid_table_alignments:
            with open(fn) as in_table:
                table_stringio = StringIO()
                table_stringio.write(in_table.read())
                table_stringio.seek(0)
                try:
                    sequence_list = dr.read_table_alignment_flo(table_stringio)
                    assert len(sequence_list) > 0
                    for i in sequence_list:
                        assert isinstance(i, Sequence)
                        assert isinstance(i.name, string_formats)
                        assert len(i.cont_values) > 0
                except Exception as e:
                    print('Raised exception: {}'.format(str(e)))
                    assert False
