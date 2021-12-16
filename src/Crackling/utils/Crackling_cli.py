from crackling import Crackling
from crackling import ConfigManager
from pathlib import Path
from crackling.Helpers import printer

import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config', help='The config file for Crackling', default=None, required=True)

    args = parser.parse_args()

    cm = ConfigManager(Path(args.config), lambda x : print(f'configMngr says: {x}'))
    if not cm.isConfigured():
        print('Something went wrong with reading the configuration.')
        exit()
    else:
        printer('Crackling is starting...')

    Crackling(cm)
    
if __name__ == '__main__':
    main()