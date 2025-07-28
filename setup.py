import os
import sys
from glob import glob
import re
import tarfile
import subprocess
import shutil
import zipfile
import tempfile
import re
import time


VERSION = "1.0.0"
AUTHOR = "PARATHAA Development Team"
MAINTAINER = "Jacob Nearing"
MAINTAINER_EMAIL = "jnearing@hsph.harvard.edu" 


# check for either of the required versions
# required python versions (3.6+)
required_python_version_major = [3]
required_python_version_minor = [6]

pass_check = False
try:
	for major, minor in zip(required_python_version_major, required_python_version_minor):
		if (sys.version_info[0] == major and sys.version_info[1] >= minor):
			pass_check = True
except (AttributeError, IndexError):
	sys.exit("CRITICAL ERROR: The python version found (version 1) " +
	         "does not match the version required (version " +
	         str(required_python_version_major) + "." +
	         str(required_python_version_minor) + "+)")

if not pass_check:
	sys.exit("CRITICAL ERROR: The python version found (version " +
	         str(sys.version_info[0]) + "." + str(sys.version_info[1]) + ") " +
	         "does not match the version required (version " +
	         str(required_python_version_major) + "." +
	         str(required_python_version_minor) + "+)")

if not pass_check:
	sys.exit("CRITICAL ERROR: The python version found (version " +
	         str(sys.version_info[0]) + "." + str(sys.version_info[1]) + ") " +
	         "does not match the version required (version " +
	         str(required_python_version_major) + "." +
	         str(required_python_version_minor) + "+)")


try:
	import setuptools
except ImportError:
	sys.exit("Please install setuptools.")
	sys.exit(1)

# check setuptools version
required_setuptools_version_major = 1
try:
	setuptools_version = setuptools.__version__
	setuptools_version_major = int(setuptools_version.split(".")[0])
	if setuptools_version_major < required_setuptools_version_major:
		sys.exit("CRITICAL ERROR: The setuptools version found (version " +
		         setuptools_version + ") does not match the version required " +
		         "(version " + str(required_setuptools_version_major) + "+)." +
					" Please upgrade your setuptools version.")
except (ValueError, IndexError, NameError):
	sys.exit("CRITICAL ERROR: Unable to call setuptools version. Please upgrade setuptools.")

from setuptools.command.install import install as _install

import distutils

# try to import urllib.request.urlretrieve for python3
try:
	from urllib.request import urlretrieve
except ImportError:
	from urllib import urlretrieve

def byte_to_megabyte(byte):
	"""
	Convert byte value to megabyte
	"""

	return byte / (1024.0 ** 2)


class ReportHook():
	def __init__(self):
		self.start_time = time.time()

	def report(self, blocknum, block_size, total_size):
		"""
		Print download progress message
		"""

		if blocknum == 0:
			self.start_time = time.time()
			if total_size > 0:
				print("Downloading file of size: " + "{:.2f}".format(byte_to_megabyte(total_size)) + " MB\n")
		else:
			total_downloaded = blocknum * block_size
			status = "{:3.2f} MB ".format(byte_to_megabyte(total_downloaded))

			if total_size > 0:
				percent_downloaded = total_downloaded * 100.0 / total_size
				# use carriage return plus sys.stdout to overwrite stdout
				try:
					download_rate = total_downloaded / (time.time() - self.start_time)
					estimated_time = (total_size - total_downloaded) / download_rate
				except ZeroDivisionError:
					download_rate = 0
					estimated_time = 0
				estimated_minutes = int(estimated_time / 60.0)
				estimated_seconds = estimated_time - estimated_minutes * 60.0
				status += "{:3.2f}".format(percent_downloaded) + " %  " + \
				          "{:5.2f}".format(byte_to_megabyte(download_rate)) + " MB/sec " + \
				          "{:2.0f}".format(estimated_minutes) + " min " + \
				          "{:2.0f}".format(estimated_seconds) + " sec "
			status += "        \r"
			sys.stdout.write(status)


def download(url, download_file):
	"""
	Download a file from a url
	"""

	try:
		print("Downloading " + url)
		file, headers = urlretrieve(url, download_file, reporthook=ReportHook().report)
		# print final return to start new line of stdout
		print("\n")
	except EnvironmentError:
		print("WARNING: Unable to download " + url)


def download_unpack_tar(url, download_file_name, folder, software_name):
	"""
	Download the url to the file and decompress into the folder
	"""

	# Check for write permission to the target folder
	if not os.access(folder, os.W_OK):
		print("WARNING: The directory is not writeable: " +
		      folder + " . Please modify the permissions.")

	download_file = os.path.join(folder, download_file_name)

	download(url, download_file)

	error_during_extract = False

	try:
		tarfile_handle = tarfile.open(download_file)
		tarfile_handle.extractall(path=folder)
		tarfile_handle.close()
	except (EnvironmentError, tarfile.ReadError):
		print("WARNING: Unable to extract " + software_name + ".")
		error_during_extract = True

	if not error_during_extract:
		try:
			os.unlink(download_file)
		except EnvironmentError:
			print("WARNING: Unable to remove the temp download: " + download_file)


def download_unpack_zip(url, download_file_name, folder, software_name):
	"""
	Download the url to the file and decompress into the folder
	"""

	# Check for write permission to the target folder
	if not os.access(folder, os.W_OK):
		print("WARNING: The directory is not writeable: " +
		      folder + " . Please modify the permissions.")

	download_file = os.path.join(folder, download_file_name)

	download(url, download_file)

	error_during_extract = False

	try:
		zipfile_handle = zipfile.ZipFile(download_file)
		zipfile_handle.extractall(path=folder)
		zipfile_handle.close()
	except EnvironmentError:
		print("WARNING: Unable to extract " + software_name + ".")
		error_during_extract = True

	if not error_during_extract:
		try:
			os.unlink(download_file)
		except EnvironmentError:
			print("WARNING: Unable to remove the temp download: " + download_file)

def download_wget_unpack_zip(url, download_file_name, folder, software_name):
	"""
	Download the url to the file and decompress into the folder
	"""

	# Check for write permission to the target folder
	if not os.access(folder, os.W_OK):
		print("WARNING: The directory is not writeable: " +
		      folder + " . Please modify the permissions.")

	download_file = os.path.join(folder, download_file_name)

	try:
		subprocess.call(["wget",url,"--directory-prefix="+folder])
	except (EnvironmentError,subprocess.CalledProcessError):				
		print("WARNING: Errors downloading mothur.")

	error_during_extract = False

	try:
		zipfile_handle = zipfile.ZipFile(download_file)
		zipfile_handle.extractall(path=folder)
		zipfile_handle.close()
	except EnvironmentError:
		print("WARNING: Unable to extract " + software_name + ".")
		error_during_extract = True

	if not error_during_extract:
		try:
			os.unlink(download_file)
		except EnvironmentError:
			print("WARNING: Unable to remove the temp download: " + download_file)


def find_exe_in_path(exe):
    """
    Check that an executable exists in $PATH
    """
    
    paths = os.environ["PATH"].split(os.pathsep)
    for path in paths:
        fullexe = os.path.join(path,exe)
        if os.path.exists(fullexe):
            if os.access(fullexe,os.X_OK):
                return path
    return None

def install_fasttree(mac_os,install_scripts, replace_install=None):
    """ 
    Download and install the FastTree software if not already installed
    """
    
    # Download the FastTree software
    # Check to see if already downloaded
    fasttree_installed=find_exe_in_path("FastTree")
    fastTree_file="FastTree"            
    fastTree_url ="http://www.microbesonline.org/fasttree/FastTree"
    download_file=os.path.join(install_scripts, fastTree_file)
    
	# install fastTree if not already installed
    if not fasttree_installed or replace_install:
        if mac_os:
            fastTree_file="FastTree.c"
            fastTree_url="http://www.microbesonline.org/fasttree/FastTree.c"
            download_file=os.path.join(install_scripts, fastTree_file)
            download(fastTree_url,download_file)
			
			#Check for GCC Version
            try:
                subprocess_output=subprocess.check_output(["gcc","--version"])
            except (EnvironmentError,subprocess.CalledProcessError):
                print("WARNING: Please install gcc.")

			#Install FastTree for MacOS	
            try:
                subprocess.call(["gcc","-O3","-finline-functions","-funroll-loops","-Wall","-o",install_scripts+"/FastTree",download_file,"-lm"])
            except (EnvironmentError,subprocess.CalledProcessError):				
                print("WARNING: Errors installing FastTree.")
        else: 
            download(fastTree_url,download_file)
            os.chmod(download_file, 0o755)
        print("Installing FastTree.")
    else:
        print("Found FastTree install.")

def install_mothur(mac_os,install_scripts, replace_install=None):
    """ 
    Download and install the Mothur software if not already installed
    """
    # Download the Mothur software
    mothur_installed=find_exe_in_path("mothur")      
    mothur_exe = "mothur"
    mothur_file="Mothur.linux_8.zip"
    mothur_url ="https://github.com/mothur/mothur/releases/download/v1.48.1/Mothur.linux_8.zip"
    
	# install fastTree if not already installed
    if not mothur_installed or replace_install:
        if mac_os:
            mothur_url="https://github.com/mothur/mothur/releases/download/v1.48.1/Mothur.OSX-10.14.zip"
            mothur_file="Mothur.OSX-10.14.zip"

        parathaa_source_folder=os.path.dirname(os.path.abspath(__file__))
        tempfolder=tempfile.mkdtemp(prefix="mothur_download_",dir=parathaa_source_folder)	
        download_wget_unpack_zip(mothur_url,mothur_file, tempfolder,mothur_exe)
        mothur_exe_full_path=os.path.join(tempfolder,mothur_exe, mothur_exe) 
               
        # copy the installed software to the final bin location
        try:
            shutil.copy(mothur_exe_full_path, install_scripts)
            os.chmod(os.path.join(install_scripts,mothur_exe), 0o755)
        except (EnvironmentError, shutil.Error):
            print("Error during Installing mothur. /bin Permission denied")
        try:
            shutil.rmtree(tempfolder)
            print("Installed mothur.")
        except EnvironmentError:
            print("WARNING: Unable to remove temp install folder.")
    else:
        print("Found mothur install.")

def install_pplacer(mac_os,install_scripts, replace_install=None):
    """ 
    Download and install the Mothur software if not already installed
    """
    # Download the pplacer software
    pplacer_installed=find_exe_in_path("pplacer")      
    pplacer_exe = "pplacer"
    pplacer_file="pplacer-linux-v1.1.alpha19.zip"
    pplacer_filename="pplacer-Linux-v1.1.alpha19"
    pplacer_url ="https://github.com/matsen/pplacer/releases/download/v1.1.alpha19/pplacer-linux-v1.1.alpha19.zip"
    
	# install pplacer if not already installed
    if not pplacer_installed or replace_install:
        if mac_os:
            pplacer_url="https://github.com/smirarab/sepp/archive/refs/tags/4.5.1.zip"
            pplacer_file="4.5.1.zip"
            pplacer_filename="sepp-4.5.1/tools/bundled/Darwin"

        parathaa_source_folder=os.path.dirname(os.path.abspath(__file__))
        tempfolder=tempfile.mkdtemp(prefix="pplacer_download_",dir=parathaa_source_folder)	
        download_wget_unpack_zip(pplacer_url,pplacer_file, tempfolder,pplacer_exe)
        pplacer_exe_full_path=os.path.join(tempfolder,pplacer_filename, pplacer_exe) 
               
        # copy the installed software to the final bin location
        try:
            shutil.copy(pplacer_exe_full_path, install_scripts)
            os.chmod(os.path.join(install_scripts,pplacer_exe), 0o755)
        except (EnvironmentError, shutil.Error):
            print("Error during Installing pplacer. /bin Permission denied")
        try:
            shutil.rmtree(tempfolder)
            print("Installed pplacer.")
        except EnvironmentError:
            print("WARNING: Unable to remove temp install folder.")
    else:
        print("Found pplacer install.")

class Install(_install):
    """
    Custom setuptools install command, set executable permissions for glpk
    """
    def run(self):
        _install.run(self)
        install_scripts=os.path.join(self.install_base, "bin")
	    # find out the platform
        mac_os=False
        if sys.platform in ["darwin","os2","os2emx"]:
            mac_os=True
        install_fasttree(mac_os,install_scripts)
        install_mothur(mac_os,install_scripts)
        install_pplacer(mac_os,install_scripts)


setuptools.setup(
	name="parathaa",
	author=AUTHOR,
	author_email=MAINTAINER_EMAIL,
	version=VERSION,
	license="MIT",
	description="PARATHAA: Preserving and Assimilating Region-specific Ambiguities in Taxonomic Hierarchical Assignments for Amplicons.",
	url="https://github.com/biobakery/parathaa",
	keywords=['microbial', 'microbiome', 'bioinformatics', 'microbiology', 'metagenomic', 'baqlava' ,'anadama2'],
	platforms=['Linux', 'MacOS'],
	classifiers=[
		"Programming Language :: Python",
		"Development Status :: 5 - Production/Stable",
		"Environment :: Console",
		"Operating System :: MacOS",
		"Operating System :: Unix",
		"Programming Language :: Python :: 3.6",
		"Topic :: Scientific/Engineering :: Bio-Informatics"
	],
	install_requires=['anadama2>=0.7.4','taxtastic'],
	packages=setuptools.find_packages(),
	cmdclass={'install': Install},
	entry_points={
		'console_scripts': [
			'parathaa_run_tree_analysis = parathaa.run_tree_analysis:main',
			'parathaa_run_taxa_assignment = parathaa.run_taxa_assignment:main',
		]},
	package_data={
		'parathaa': [
			'utility/*',
			'data/*',
		]},
	scripts=['parathaa/parathaa_plot_assignment.R'],
	zip_safe=False
)
