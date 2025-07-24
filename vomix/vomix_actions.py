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
        # sf = os.path.join(os.path.dirname(os.path.realpath("vomix/workflow/rules/")), filename)
        sf = os.path.realpath("vomix/workflow/rules/" + filename)
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
            if value is not None and attr != 'custom_config' and attr != 'name':
                attr = str.replace(attr, "_", "-")
                script += f'{attr}="{value}" '

        for attr, value in snakemake_obj.__dict__.items():
            if value is not None and attr != 'add_args':
                attr = str.replace(attr, "_", "-")
                script += f'{attr}="{value}" '

        script += "--sdm conda --use-conda"

        return script

    def createFoldersAndUpdateConfig(self, module_obj):
        # Create outdir + datadir folders 
        outdir = module_obj.outdir
        datadir = module_obj.datadir

        if not (os.path.exists(outdir) and os.path.exists(os.path.join(outdir, ".vomix"))):
            os.makedirs(os.path.join(outdir, ".vomix"), exist_ok=True)

        now = datetime.datetime.now()
        latest_run = now.strftime("%Y%m%d_%H%M%S")
        outdir_folder = os.path.join(outdir, ".vomix/log/vomix" + latest_run)
        datadir_folder = datadir

        os.makedirs(outdir_folder, exist_ok=True)
        os.makedirs(datadir_folder, exist_ok=True)

        # if custom config is specified
        if module_obj.custom_config is not None:
            logging.info(f"Using custom config: {module_obj.custom_config}")
            logging.info(f"REMINDER - Any command line flags spefied will override those options in your custom config.")
            shutil.copy(os.path.realpath(module_obj.custom_config), outdir_folder)
            os.rename(outdir_folder + "/" + module_obj.custom_config, outdir_folder + "/config.yml")
        else:
            # Create a new config file from the config template
            logging.info(f"Using template config: config/config.yml")
            shutil.copy(os.path.realpath("config/config.yml"), outdir_folder)
    
        # edit new config with user options + latest_run
        with open(outdir_folder + "/config.yml") as f:
            list_doc = yaml.safe_load(f)
            list_doc["latest-run"] = latest_run

            for module in module_obj.__dict__:
                value = module_obj.__dict__[module]
                if value is not None and module != 'custom_config':
                    module = str.replace(module, "_", "-")
                    list_doc[module] = value
                    # logging.info(f"///////TEST: {module} // " , {value})


        with open(outdir_folder + "/config.yml", "w") as f:
            yaml.dump(list_doc, f)

        return outdir_folder

    def run_module(self, module, module_obj, snakemake_obj):

        outdir_folder = self.createFoldersAndUpdateConfig(module_obj)

        # create the script to run the module 
        script_path = os.path.realpath(outdir_folder + "/snakemake" +".sh")
        script = self.createScript(module, module_obj, snakemake_obj)
        
        # save the script
        # TODO run .sh in vomix
        script_path = os.path.realpath("snakemake" +".sh")
        # script_path = str.replace(abspath(getsourcefile(lambda:0)), "vomix_actions.py", "snakemake.sh")
        # logging.info(f"///////////////// script: {script_path}")
        with open(script_path, "w") as f:
            f.write(script)

        # Run script
        logging.info(f"Running script: {script_path}")
        cmd = ['bash', script_path]

        try:
            with Popen(cmd, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
                for line in p.stdout:
                    print(line, end='') 
            if p.returncode != 0:
                raise CalledProcessError(p.returncode, p.args)

        except subprocess.CalledProcessError as e:
            return f"Error: {e.stderr}"
        
    # def run_last_module(module) -> str:
    #     try:
    #         script_path = os.path.realpath("vomix/snakemake.sh")

    #         with open(script_path, "r") as file:
    #             lines = file.readlines()

    #         if not lines:
    #             return f"Error: The script file is empty."

    #         with open(script_path, "w") as file:
    #             for line in lines:
    #                 if 'snakemake --config module=' in line:
    #                     # Replace the module value in the line
    #                     parts = line.split('module="')
    #                     if len(parts) > 1:
    #                         rest = parts[1].split('"', 1)
    #                         if len(rest) > 1:
    #                             line = parts[0] + f'module="{module}' + '"' + rest[1][rest[1].find('"'):]
    #                 file.write(line)

    #         cmd = ['bash', script_path]

    #         # result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    #         result = subprocess.run(cmd, stdout=open('out.log', 'w'), stderr=open('error.log', 'a'),)
    #         return result.stdout 
    #     except subprocess.CalledProcessError as e:
    #         return f"Error: {e.stderr}"
