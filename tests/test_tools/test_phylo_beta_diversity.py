"""This module is used for testing phylogenetic beta diversity."""
import os

import numpy as np

from lmpy import Matrix

from biotaphy.tools.phylo_beta_diversity import cli


# .....................................................................................
def _check_csv_matrix_equality(
    csv_filename_1,
    csv_filename_2,
    num_header_rows=1,
    num_header_cols=1,
):
    """Check if two csv matrix files are approximately equal.

    Args:
        csv_filename_1 (str): Path to the first csv matrix file.
        csv_filename_2 (str): Path to the second csv matrix file.
        num_header_rows (int): The number of header rows in the matrices.
        num_header_cols (int): The number of header columns in the matrices.

    Returns:
        bool: Indication if the two csv matrices are approximately equal.
    """
    with open(csv_filename_1, mode='rt') as in_csv_1:
        mtx_1 = Matrix.load_csv(
            in_csv_1, num_header_rows=num_header_rows, num_header_cols=num_header_cols
        )
    with open(csv_filename_2, mode='rt') as in_csv_2:
        mtx_2 = Matrix.load_csv(
            in_csv_2, num_header_rows=num_header_rows, num_header_cols=num_header_cols
        )
    return np.allclose(mtx_1, mtx_2)


# .....................................................................................
def test_phylo_beta_diversity_jaccard(
    monkeypatch,
    tmpdir,
    valid_phylo_beta_diversity_package,
):
    """Test phylogenetic beta diversity using the jaccard index family.

    Args:
        monkeypatch (pytest.fixture): A fixture for monkeypatching the environment.
        tmpdir (pytest.fixture): A fixture that provides a temporary directory.
        valid_phylo_beta_diversity_package (tuple): A tuple of information that
            together forms a valid phylogenetic beta diversity package.

    Note:
        Test values were determined from example at
            https://rdrr.io/rforge/betapart/man/phylo.beta.pair.html
    """
    (pam_fn, tree_fn, test_beta_jac_fn, test_beta_jne_fn, test_beta_jtu_fn,
     _, _, _, test_phylo_beta_jac_fn, test_phylo_beta_jne_fn,
     test_phylo_beta_jtu_fn, _, _, _) = valid_phylo_beta_diversity_package

    data_format = 'csv'
    out_dir = tmpdir.strpath

    params = [
        'phylo_beta_diversity.py',
        tree_fn,
        pam_fn,
        data_format,
        'jaccard',
        out_dir,
    ]

    monkeypatch.setattr('sys.argv', params)
    cli()

    assert _check_csv_matrix_equality(
        test_beta_jac_fn, os.path.join(out_dir, 'beta_jac.csv')
    )
    assert _check_csv_matrix_equality(
        test_beta_jne_fn, os.path.join(out_dir, 'beta_jne.csv')
    )
    assert _check_csv_matrix_equality(
        test_beta_jtu_fn, os.path.join(out_dir, 'beta_jtu.csv')
    )
    assert _check_csv_matrix_equality(
        test_phylo_beta_jac_fn, os.path.join(out_dir, 'phylo_beta_jac.csv')
    )
    assert _check_csv_matrix_equality(
        test_phylo_beta_jne_fn, os.path.join(out_dir, 'phylo_beta_jne.csv')
    )
    assert _check_csv_matrix_equality(
        test_phylo_beta_jtu_fn, os.path.join(out_dir, 'phylo_beta_jtu.csv')
    )


# .....................................................................................
def test_phylo_beta_diversity_sorensen(
    monkeypatch,
    tmpdir,
    valid_phylo_beta_diversity_package,
):
    """Test phylogenetic beta diversity using the sorensen index family.

    Args:
        monkeypatch (pytest.fixture): A fixture for monkeypatching the environment.
        tmpdir (pytest.fixture): A fixture that provides a temporary directory.
        valid_phylo_beta_diversity_package (tuple): A tuple of information that
            together forms a valid phylogenetic beta diversity package.

    Note:
        Test values were determined from example at
            https://rdrr.io/rforge/betapart/man/phylo.beta.pair.html
    """
    (pam_fn, tree_fn, _, _, _, test_beta_sim_fn, test_beta_sne_fn,
     test_beta_sor_fn, _, _, _, test_phylo_beta_sim_fn,
     test_phylo_beta_sne_fn, test_phylo_beta_sor_fn
     ) = valid_phylo_beta_diversity_package

    data_format = 'csv'
    out_dir = tmpdir.strpath

    params = [
        'phylo_beta_diversity.py',
        tree_fn,
        pam_fn,
        data_format,
        'sorensen',
        out_dir,
    ]
    print(params)
    monkeypatch.setattr('sys.argv', params)
    cli()

    assert _check_csv_matrix_equality(
        test_beta_sim_fn, os.path.join(out_dir, 'beta_sim.csv')
    )
    assert _check_csv_matrix_equality(
        test_beta_sne_fn, os.path.join(out_dir, 'beta_sne.csv')
    )
    assert _check_csv_matrix_equality(
        test_beta_sor_fn, os.path.join(out_dir, 'beta_sor.csv')
    )
    assert _check_csv_matrix_equality(
        test_phylo_beta_sim_fn, os.path.join(out_dir, 'phylo_beta_sim.csv')
    )
    assert _check_csv_matrix_equality(
        test_phylo_beta_sne_fn, os.path.join(out_dir, 'phylo_beta_sne.csv')
    )
    assert _check_csv_matrix_equality(
        test_phylo_beta_sor_fn, os.path.join(out_dir, 'phylo_beta_sor.csv')
    )
