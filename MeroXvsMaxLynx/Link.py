## This class keeps crosslinking info within dictionaries that contains protein:key its link-site:values
## There are two dictionaries for di-peptide crosslinking
import json


class Link:


    def __init__(self, xlsites1, xlsites2):
        self.xlsites1 = xlsites1
        self.xlsites2 = xlsites2


    def __eq__(self, other):
        if isinstance(other, Link):
            return self.xlsites1 == other.xlsites1 and self.xlsites2 == other.xlsites2
        return False

    def __str__(self):
        return json.dumps(self.xlsites1) + json.dumps(self.xlsites2)
