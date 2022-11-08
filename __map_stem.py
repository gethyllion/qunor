from  argparse import ArgumentParser 
from qunor import normalize

Parser = ArgumentParser(description = normalize.stemness_note)
Parser.add_argument('-n','--normalize',type=str)

cmd_parser = Parser.parse_args()

df = cmd_parser.normalize
m = normalize.Normalize(df, 
                        only_stemness=True)
m.stemness(save=True)
