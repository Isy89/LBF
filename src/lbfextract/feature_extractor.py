from lbfextract.pluggin_manager import get_pluggin_manager


def flatten_list(l):
    """function to flatten a list in a recursive way"""
    if isinstance(l, list):
        return [a for i in l for a in flatten_list(i)]
    else:
        return [l]


class FeatureExtractor:
    """
    This object provide a python interface to all feature extraction methods implemented. 
    All commands can be used from python trhough this object.
    """

    def __init__(self):
        pm = get_pluggin_manager()
        self.extractors = {i.name.replace("-", "_"): i for i in flatten_list(pm.hook.get_command())}

    def extract(self, extractor_name: str, **kwargs):
        dict_params = {param.name: param.default for param in self.extractors[extractor_name].params}
        dict_params.update(kwargs)
        return self.extractors[extractor_name].callback(**dict_params)

    def get_help_for_extractor(self, extractor_name):
        backslah_n = "\n"
        params = ' '.join([
            f"    {param.name}({str(param.default)}) => {str(param.help)} {backslah_n}"
            for param in self.extractors[extractor_name].params
        ])
        msg = f"extractor {extractor_name} with following parameters:\n{params}"
        print(msg)

    def get_exctractor_names(self):
        return list(self.extractors.keys())

    def help(self):
        msg = []
        backslah_n = "\n"
        for i in self.extractors:
            params = ' '.join([
                f"    {param.name}({str(param.default)}) => {str(param.help)}  {backslah_n}"
                for param in self.extractors[i].params
            ])
            msg.append(
                f"extractor {i} with following parameters:\n{params}"
            )
        return '\n'.join(msg)

    def __str__(self):
        return self.help()


feature_extractor = FeatureExtractor()

