from pathlib import Path
import argparse

from haplocrackling import HaploCrackling
from haplocrackling import ConfigManager
from haplocrackling.Helpers import printer

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config', help='The config file for HaploCrackling', default=None, required=True)

    args = parser.parse_args()

    cm = ConfigManager(Path(args.config), lambda x : print(f'configMngr says: {x}'))
    if not cm.isConfigured():
        print('Something went wrong with reading the configuration.')
        exit()
    else:
        printer('HaploCrackling is starting...')

    HaploCrackling(cm)
    
if __name__ == '__main__':
    main()