
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri

class AcousticIndices:
  def __init__(self, script_path = "indices.R"):
    self.__R = robjects.r
    self.__R.source(script_path)
    pandas2ri.activate() #automatic convertion from R dataframe to Pandas dataframe 

  def process_file(self, path, aquatic = False, slice_size = 1):
    return self.__R["process.file"](path = path,
                                    aquatic = ("TRUE" if aquatic else "FALSE"),
                                    slice_size = slice_size,
                                    start_parallel = "TRUE")

  def process_dir(self, source_path, target_path = None, aquatic = False, slice_size = 1):
    return self.__R["process.dir"](source_path = source_path,
                                   target_path = source_path if target_path is None else target_path,
                                   aquatic = ("TRUE" if aquatic else "FALSE"),
                                   slice_size = slice_size)
                             
