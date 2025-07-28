import subprocess
import os
import sys
from subprocess import Popen, PIPE, CalledProcessError
import datetime
import json
import yaml
import shutil
import logging
from inspect import getsourcefile
from os.path import abspath
import inspect, os.path
import re


logging.basicConfig(level=logging.INFO)

class vomix_actions:
    def __init__(self):
        self.name = "vomix"
        self.version = "1.0.0"
        self.description = "vomix is a tool for viral metagenomics analysis."

    def __repr__(self):
        return f"vomix(name={self.name}, version={self.version}, description={self.description})"

    def __str__(self):
        return f"{self.name} v{self.version}: {self.description}"
    
    def get_snakefile(filename):
        filename = inspect.getframeinfo(inspect.currentframe()).filename
        sf     = str.replace(os.path.dirname(os.path.abspath(filename)), "/.venv/lib/python3.9/site-packages/vomix", "vomix/workflow/Snakefile")

        if not os.path.exists(sf):
            sys.exit("Unable to locate the Snakemake file; tried %s" % sf)
        return sf
    
    def env_setup_script() -> str:

        script_path = os.path.realpath("vomix/env_setup.sh")
        logging.info(f"Running script: {script_path}")

        cmd = ['bash', script_path]
        try:
            result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            return result.stdout
        except subprocess.CalledProcessError as e:
            return f"Error: {e.stderr}"
    

    def createScript(self, module, module_obj, snakemake_obj):
        script = ""

        script += "snakemake --config module=\"" + module + "\" "

        for attr, value in module_obj.__dict__.items():
            attr = str.replace(attr, "_", "-")
            if str(value) == 'True' or str(value) == 'False':
                script += f'{attr}={value} '
            elif value is not None and attr != 'custom_config' and attr != 'name':
                attr = str.replace(attr, "_", "-")
                script += f'{attr}="{value}" '


        for attr, value in snakemake_obj.__dict__.items():
            if value is not None and attr != 'add_args':
                attr = str.replace(attr, "_", "-")
                if str(value) == 'True':
                    script += f'--{attr} '
                elif str(attr) == 'jobs':
                    script += f'--{attr} {value} '
                elif str(attr) == 'latency-wait':
                    script += f'--{attr} {value} '
                elif str(value) != 'False':
                    script += f'{attr}="{value}" '
            if attr == 'add_args' and value is not None and value != '':
                script += f'{value} '

        script += "--sdm conda --use-conda"

        return script

    def createFoldersAndUpdateConfig(self, module_obj):
        currentVomixDir = os.path.dirname(os.path.abspath(__file__))
        count = sum(1 for _ in re.finditer(r'\b%s\b' % re.escape("/vomix"), currentVomixDir))
        configPath = "/config/config.yml"
        fullConfigPath = ""

        if count == 1:
            fullConfigPath = str.replace(currentVomixDir, "/vomix", configPath)
        elif count > 1:
            fullConfigPath = configPath.join(currentVomixDir.rsplit("/vomix", count - 1))
        else:
            raise FileNotFoundError("Could not determine the path to the template config file.")
        
        logging.info(f"Using template config: {fullConfigPath}")

        # get workdir
        workdir = module_obj.workdir
        if workdir is None:
            # get workdir from custom config
            if module_obj.custom_config is not None:
                with open(module_obj.custom_config) as f:
                    list_doc = yaml.safe_load(f)
                    workdir = list_doc["workdir"]
            else:
                # get workdir from template config.yml 
                with open(fullConfigPath) as f:
                    list_doc = yaml.safe_load(f)
                    workdir = list_doc["workdir"]
        
        logging.info(f"Using workdir: {workdir}")

        # Create outdir + datadir folders 
        outdir = workdir + module_obj.outdir
        if module_obj.datadir is not None: 
            datadir = workdir + module_obj.datadir
            datadir_folder = datadir
            os.makedirs(datadir_folder, exist_ok=True)
        if module_obj.fastadir is not None: 
            fastadir = workdir + module_obj.fastadir
            fastadir_folder = fastadir
            os.makedirs(fastadir_folder, exist_ok=True)

        if not (os.path.exists(outdir) and os.path.exists(os.path.join(outdir, ".vomix"))):
            os.makedirs(os.path.join(outdir, ".vomix"), exist_ok=True)

        now = datetime.datetime.now()
        latest_run = now.strftime("%Y%m%d_%H%M%S")
        outdir_folder = os.path.join(outdir, ".vomix/log/vomix" + latest_run)

        os.makedirs(outdir_folder, exist_ok=True)

        # if datadir_folder is not None and datadir_folder != "":
        #     os.makedirs(datadir_folder, exist_ok=True)
        # if fastadir_folder is not None and fastadir_folder != "":
        #     os.makedirs(fastadir_folder, exist_ok=True)
            
        # if custom config is specified
        if module_obj.custom_config is not None:
            logging.info(f"Using custom config: {module_obj.custom_config}")
            logging.info(f"REMINDER - Any command line flags spefied will override those options in your custom config.")
            shutil.copy(os.path.realpath(module_obj.custom_config), outdir_folder)
            os.rename(outdir_folder + "/" + module_obj.custom_config, outdir_folder + "/config.yml")
        else:
            # Create a new config file from the config template

            # currentVomixDir = os.path.dirname(os.path.abspath(__file__))
            # count = sum(1 for _ in re.finditer(r'\b%s\b' % re.escape("/vomix"), currentVomixDir))
            # configPath = "/config/config.yml"

            # if count == 1:
            #     path = str.replace(currentVomixDir, "/vomix", configPath)
            # elif count > 1:
            #     path = configPath.join(currentVomixDir.rsplit("/vomix", count - 1))
            # else:
            #     raise FileNotFoundError("Could not determine the path to the template config file.")

            # logging.info(f"Using template config: {path}")
    
            shutil.copy(fullConfigPath, outdir_folder)

        # edit new config with user options + latest_run
        with open(outdir_folder + "/config.yml") as f:
            list_doc = yaml.safe_load(f)
            list_doc["latest-run"] = latest_run

            for module in module_obj.__dict__:
                value = module_obj.__dict__[module]
                if value is not None:
                    module = str.replace(module, "_", "-")
                    list_doc[module] = value
                    logging.debug(f"Setting config option: {module} = {value}")
               
        with open(outdir_folder + "/config.yml", "w") as f:
            yaml.dump(list_doc, f)

        return outdir_folder

    def run_module(self, module, module_obj, snakemake_obj):

        outdir_folder = self.createFoldersAndUpdateConfig(module_obj)

        # create the script to run the module 
        script_path = os.path.realpath(outdir_folder + "/snakemake" +".sh")
        script = self.createScript(module, module_obj, snakemake_obj)
        
        # save the script
        with open(script_path, "w") as f:
            f.write(script)

        logging.info(f"Running script: {script_path}")
        cmd = ['bash', script_path]

        currentVomixDir = os.path.dirname(os.path.abspath(__file__))
        count = sum(1 for _ in re.finditer(r'\b%s\b' % re.escape("/vomix"), currentVomixDir))
        upPath = "/"

        if count == 1:
            currentWorkingPath = str.replace(currentVomixDir, "/vomix", upPath)
        elif count > 1:
            currentWorkingPath = upPath.join(currentVomixDir.rsplit("/vomix", count - 1))
        else:
            raise FileNotFoundError("Could not determine the path to the current working directory.")

        logging.info(f"currentWorkingPath: {currentWorkingPath}")

        try:
            with Popen(cmd, stdout=PIPE, bufsize=1, universal_newlines=True, cwd=currentWorkingPath) as p:
                for line in p.stdout:
                    print(line, end='') 

            if p.returncode != 0:
                raise CalledProcessError(p.returncode, p.args)


        except subprocess.CalledProcessError as e:
            return f"Error: {e.stderr}"
        