import re
import pickle
import pathlib
from sklearn.model_selection import train_test_split


class GenerateTypes(object):
    """Generates types files for the training and validation sets for EGNN regression.
    
    Noes:
        EGNN requires 2 major pre-preperation steps; generating types files and generating parquet files.
        The types files contain assay data (labels) and paths to the parquet files.

    Args:
        index_file (file): pdbbind index file
        id_file (file): a txt file containing the pdb ids for the desired proteins (1 id per line)
        random_seed (int): the random seed index for train/validation splitting
        mode (bool): if not False, must be specified as either 'redocking' or 'crossdocking'
        crossdock_dict_file: ???
        file_string: ???
        out_path (path): relative path to the directory in which the types files should be stored

    Returns:
        Generates 2 files; a training set types file and a testing set types file

    Examples:
        >>> GenerateTypes(
            '../Data/pdbbind_2020_general/index/INDEX_general_PL_data.2020', 
            '../Data/ids/filtered_pdbbind_2020_general_id.txt',
            out_path='../Data/types/'
            )
    
    """
    def __init__(
        self,
        index_file,
        id_file,
        random_seed=42,
        mode=False,
        crossdock_dict_file=None,
        file_string=None,
        out_path=None,
    ):

        #check files exist
        assert pathlib.Path(index_file).is_file(), "File does not exist"
        assert pathlib.Path(id_file).is_file(), "File does not exist"
        if crossdock_dict_file is not None:
            assert pathlib.Path(crossdock_dict_file).is_file(), "File does not exist"
        if out_path is not None:
            assert pathlib.Path(out_path).is_dir(), "Directory does not exist"
        
        #check seed is an integer
        assert isinstance(random_seed, int), "Random seed must be an integer"

        #check mode is legal
        assert mode in [False, 'redocking', 'crossdocking'], "The mode is not a legal option"

        self.index_file = index_file
        self.mode = mode
        self.crossdock_dict_file = crossdock_dict_file
        self.file_string = file_string

        self.out_path = out_path

        self.pdb_list = GenerateTypes._lines_to_list(id_file)

        self.train_list, self.val_list = GenerateTypes._train_val_split(
            self.pdb_list, random_seed
        )

        self.train_str, self.val_str = self.generate_types()

        self.write_types()

    def generate_types(self):

        train_str = ""
        val_str = ""

        with open(self.index_file) as f:
            lines = f.readlines()[6:]
            for line in lines:
                pdb = line.split()[0]
                if pdb in self.pdb_list:
                    types_str = ""
                    protein_parquet_file = f"{pdb}/{pdb}_protein_cleaned.parquet"
                    if not self.mode:
                        ligand_parquet_file = f"{pdb}/{pdb}_ligand_pymol.parquet"
                    elif self.mode == "redocking":
                        ligand_parquet_file = f"{pdb}/{pdb}_redocked_best_pose.parquet"
                    elif self.mode == "crossdocking":
                        cross_dock_dict = pickle.load(
                            open(self.crossdock_dict_file, "rb")
                        )
                        ligand_parquet_file = f"{pdb}/{pdb}_crossdocking_{self.file_string}_similar_best_pose.parquet"
                        protein_parquet_file = f"{cross_dock_dict[pdb]}/{cross_dock_dict[pdb]}_protein_cleaned.parquet"
                    affinity, metric = None, None
                    affinity = line.split()[3]
                    metric = (re.split(r"~|=|>|<", line.split()[4]))[0]
                    types_line = GenerateTypes._make_types_line(
                        protein_parquet_file, ligand_parquet_file, affinity, metric
                    )
                    types_str += types_line
                    if pdb in self.train_list:
                        train_str += types_str
                    elif pdb in self.val_list:
                        val_str += types_str

        return train_str[:-1], val_str[:-1]

    def write_types(self):

        with open(self.out_path + "pdbbind_2020_general_pymol_train.types", "w") as f:
            f.write(self.train_str)
        with open(self.out_path + "pdbbind_2020_general_pymol_val.types", "w") as f:
            f.write(self.val_str)

    @staticmethod
    def _lines_to_list(file):
        if file == None:
            return []
        else:
            with open(file) as f:
                cleaned_list = [line.strip() for line in f.readlines()]
            return cleaned_list

    @staticmethod
    def _make_types_line(receptor_pdb, ligand_sdf, affinity, metric):
        affinities = [-1, -1, -1]
        affinities[["Ki", "Kd", "IC50"].index(metric)] = affinity
        affinity_str = "{0} {1} {2}".format(*affinities)
        return "{0} {1} {2}\n".format(affinity_str, receptor_pdb, ligand_sdf)

    @staticmethod
    def _train_val_split(pdb_indices, r_seed):

        pdb_array = pdb_indices
        train, val = train_test_split(pdb_array, random_state=r_seed, test_size=0.2)
        return train, val
