#!/usr/bin/env python

'''
@Copyright: Peking University Cancer Hosptial, All Rights Reserved.
@Author: Du Yang
@Date: 2020-02-28 11:03:40
@LastEditors: Du Yang
@LastEditTime: 2020-03-12 10:15:23
@Description: create, update, and combine config file
'''

import configparser
import os
import sys
import argparse
import getpass
from collections import Counter

if sys.getdefaultencoding() != 'utf-8':
    reload(sys)
    sys.setdefaultencoding('utf-8')
    
def read_config(config):
    cf = configparser.RawConfigParser()
    cf.optionxform = str
    cf.read(config)
    return(cf)

def update_config(config, section, key, value):
    
    if not os.path.exists(config):
        print("config file " + config + " is not exist")
        print("create new config file " + config + "")
        open(config, 'a').close()

    cf = read_config(config)

    if not cf.has_section(section):
        cf.add_section(section)
    
    if cf.has_option(section, key):
        print("update existed key: " + key + " in section: " + section)
        print(cf.get(section, key) + " ====> " + value)

    else:
        print("add new key: " + key + " in section: " + section)
        print(key + " ====> " + value)
    
    cf.set(section, key, value)

    with open(config, 'w') as fw:
        cf.write(fw)

def handle_args():

    usage = "create and update config file"
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('config_file', default='config.ini', help="config file", metavar="")
    parser.add_argument('section', default='Default', help="config section", metavar="")
    parser.add_argument('key', default='test', help="config key", metavar="")
    parser.add_argument('value', default='test', help="config value", metavar="")
    args = parser.parse_args(sys.argv[1:])
    return args

def main(options):

    config_file = options.config_file
    section = options.section
    key = options.key
    value = options.value

    update_config(config_file, section, key, value)

if __name__ == "__main__":
    options = handle_args()
    main(options)
