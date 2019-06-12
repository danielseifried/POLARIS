#!/usr/bin/env bash

# Set Install directory
install_directory=$(pwd)

fits_support=false
polaris_tools=false

function command_exists {
    type "$1" &> /dev/null ;
}

# ---------------------------------------------------------------------------------
# ------------------------- Routines for installations ----------------------------
# ---------------------------------------------------------------------------------
function install_fits_support {
    # Check for necessary programms
    if ! command_exists cmake
    then
        echo -e "--- ${RED}Error:${NC} cmake not found. Please install cmake!"
        exit
    fi
    # Install Libraries
    echo -e "${PC}--- Install required libraries for fits support ---${NC}"
    echo -e "Install cfitsio"
    if [ ! -d "cfitsio" ]
    then
        tar -xf cfitsio.tar.gz
    fi
    cd cfitsio
    echo -ne  "- Configuring cfitsio ... "\\r
    cmake . -DBUILD_SHARED_LIBS=OFF > /dev/null 2>&1 \
        && echo -e "- Configuring cfitsio [${GREEN}done${NC}]" \
        || { echo -e  "- Configuring cfitsio [${RED}Error${NC}]"; exit; }
    echo -ne  "- Compiling cfitsio ... "\\r
    make > /dev/null 2>&1 \
        && echo -e "- Compiling cfitsio [${GREEN}done${NC}]" \
        || { echo -e  "- Compiling cfitsio [${RED}Error${NC}]"; exit; }
    cd ..
    
    echo -e "Install CCfits"
    if [ ! -d "CCfits" ]
    then
        tar -xf CCfits.tar.gz
    fi
    cd CCfits
    echo -ne "- Configuring CCfits ... "\\r
    cmake . -DCMAKE_PREFIX_PATH=${install_directory}/cfitsio > /dev/null 2>&1 \
        && echo -e "- Configuring CCfits [${GREEN}done${NC}]" \
        || { echo -e  "- Configuring CCfits [${RED}Error${NC}]"; exit; }
    echo -ne "- Compiling CCfits ... "\\r
    make > /dev/null 2>&1 \
        && echo -e "- Compiling CCfits [${GREEN}done${NC}]" \
        || { echo -e  "- Compiling CCfits [${RED}Error${NC}]"; exit; }
    cd ..

    export_str="export LD_LIBRARY_PATH=\"${install_directory}/CCfits/.libs/:${install_directory}/cfitsio:"'${LD_LIBRARY_PATH}'"\""
    if ! grep -q "${export_str}" ${HOME}/.bashrc
    then
            echo "${export_str}" >> ${HOME}/.bashrc
            echo -e "- Updating bashrc [${GREEN}done${NC}]"
    fi
    source ~/.bashrc
}

function check_python_packages {
    echo -e  "Looking for required Python packages"
    for package_name in numpy scipy matplotlib astropy os array struct argparse
    do
        if ! python -c "import ${package_name}" &> /dev/null
        then
            echo -e  "- Required python package ${package_name} [${RED}not found${NC}]"
            if command_exists conda && [ ${package_name} != "matplotlib2tikz" ]
            then
                echo "--- Conda installation detected. Install ${package_name}!"
                conda install ${package_name}
                if ! python -c "import ${package_name}" &> /dev/null
                then
                    echo -e  "- ${RED}Error:${NC} installation of python packages not succesfull!"
                    exit
                else
                    echo -e  "- Required python package ${package_name} [${GREEN}found${NC}]"
                fi
            elif command_exists pip
            then
                echo "--- Pip installation detected. Install ${package_name}!"
                pip install ${package_name}
                if ! python -c "import ${package_name}" &> /dev/null
                then
                    echo -e  "- ${RED}Error:${NC} installation of python packages not succesfull!"
                    exit
                else
                    echo -e  "- Required python package ${package_name} [${GREEN}found${NC}]"
                fi
            else
                echo -e "--- ${RED}Error:${NC} Please install the required python packages via package manager!"
                exit
            fi
        else
            echo -e  "- Required python package ${package_name} [${GREEN}found${NC}]"
        fi
    done
}

function install_anaconda {
    unameOut="$(uname -s)"
    case "${unameOut}" in
        Linux*)  
            MACHINE_TYPE=`uname -m`
            if [ "${MACHINE_TYPE}" == 'x86_64' ]; then
                if [ ! -f "Miniconda3-latest-Linux-x86_64.sh" ]
                then
                    wget "https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh"
                fi
                chmod +x "Miniconda3-latest-Linux-x86_64.sh"
                ./Miniconda3-latest-Linux-x86_64.sh
                rm "Miniconda3-latest-Linux-x86_64.sh"
            else
                if [ ! -f "Miniconda3-latest-Linux-x86.sh" ]
                then
                    wget "https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86.sh"
                fi
                chmod +x "Miniconda3-latest-Linux-x86.sh"
                ./Miniconda3-latest-Linux-x86.sh
                rm "Miniconda3-latest-Linux-x86.sh"
            fi
        ;;
        Darwin*)
            if [ ! -f "Miniconda3-latest-MacOSX-x86_64.sh" ]
            then
                wget "https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh"
            fi
            chmod +x "Miniconda3-latest-MacOSX-x86_64.sh"
            ./Miniconda3-latest-MacOSX-x86_64.sh
            rm "Miniconda3-latest-MacOSX-x86_64.sh"
        ;;
        *)
        echo "Machine not known!"
        ;;
    esac
    source ~/.bashrc
}

function install_polaris_tools {
    # Check for installation of python
    echo -ne "Checking for Python installation ...${NC}"\\r
    python_version="$(python -V 2>&1)"
    if [ "${python_version:7:1}" -lt 3 ] || [ "${python_version:9:1}" -lt 5 ]
    then
        echo -e "${RED}Error:${NC} No python >= 3.5 installation found!"
        printf '%s\n' "Do you want to install anaconda python package? [Y/n]  "
        read install_anaconda_ans
        case ${install_anaconda_ans:=y} in
            [yY]*) 
                if ! command_exists wget
                then
                    echo -e "--- ${RED}Error:${NC} Programm wget is not installed on this system!"
                    exit
                fi   
                install_anaconda
                echo -e "Checking for Python installation [${GREEN}found${NC}]"
                check_python_packages
            ;;

            *) 
                server_support=false
            ;; 
        esac
        python_installed=false
    else
        echo -e "Checking for Python installation [${GREEN}found${NC}]"
        python_installed=true
        check_python_packages
    fi
    export_str="export PATH=\"/home/rbrauer/.local/bin:"'$PATH'"\""
    if grep -q "${export_str}" ${HOME}/.bashrc
    then
            true
    else
            echo -ne "Updating bashrc ..."\\r
            echo "${export_str}" >> ${HOME}/.bashrc
            echo -e "Updating bashrc [${GREEN}done${NC}]"
    fi
}

function install_polaris {
    # Check for necessary programms
    if ! command_exists cmake
    then
        echo -e "--- ${RED}Error:${NC} cmake not found. Please install cmake!"
        exit
    fi
    sed -i.bak 's,^//#define FITS_EXPORT,#define FITS_EXPORT,g' "${install_directory}/src/typedefs.h"
    cmake . -DBUILD_STATIC_LIBS=$1 > /dev/null 2>&1 \
        && echo -e "Configuring POLARIS [${GREEN}done${NC}]" \
        || { echo -e  "Configuring POLARIS [${RED}Error${NC}]"; exit; }
    echo -ne "Compiling POLARIS ... "\\r
    make > /dev/null 2>&1 \
        && echo -e "Compiling POLARIS [${GREEN}done${NC}]" \
        || { echo -e  "Compiling POLARIS [${RED}Error${NC}]"; exit; }
    echo -ne "Installing POLARIS ... "\\r
    make install > /dev/null 2>&1 \
        && echo -e "Installing POLARIS [${GREEN}done${NC}]" \
        || { echo -e  "Installing POLARIS [${RED}Error${NC}]"; exit; }

    export_str="export PATH=\"${install_directory}/bin:"'$PATH'"\""
    if grep -q "${export_str}" ${HOME}/.bashrc
    then
            true
    else
            echo -ne "Updating bashrc ..."\\r
            echo "${export_str}" >> ${HOME}/.bashrc
            echo -e "Updating bashrc [${GREEN}done${NC}]"
    fi
}

# Define colors
RED='\033[0;31m'
GREEN='\033[0;32m'
PC='\033[0;35m'
TC='\033[0;36m'
NC='\033[0m'

function usage {
    echo "usage: install_polaris.sh [-u] [-h]"
    echo ""
    echo "Install of POLARIS on your system!"
    echo ""
    echo "optional arguments:"
    echo "-h      show this help message and exit"
    echo "-r      clean and compile POLARIS including PolarisTools if enabled (release mode)"
    echo "-d      clean and compile POLARIS including PolarisTools if enabled (debug mode)"
    echo "-u      compile POLARIS including PolarisTools if necessary (using release mode)"
    echo "-c      pre-compile POLARIS to be portable on other machines (using release mode)"
    echo "-D      delete POLARIS from your computer"
    exit
}

while getopts "hrducD" opt; do
    case $opt in
	h)
	    usage
        ;;
	r)
	    echo -e "${TC}------ clean and compile POLARIS (${GREEN}release mode!${TC}) ------${NC}"
	    cd ${install_directory}
	    make clean
        make && make install
        if [ -d "tools" ]
        then
	        cd "tools"
	    fi
	    python setup.py install --user &> /dev/null
	    exit
	    ;;
	d)
            echo -e "${TC}------ clean and compile POLARIS (${RED}debug mode!${TC}) ------${NC}"
	    cd ${install_directory}
	    make clean
        make CXXFLAGS='-O0 -g'  && make install #'-Wall -Weffc++ -Wextra -Wsign-conversion'
        if [ -d "tools" ]
        then
	        cd "tools"
	    fi
	    python setup.py install --user &> /dev/null
	    exit
	    ;;
    u)
	    echo -e "${TC}------ compile POLARIS ------${NC}"
        cd ${install_directory}
        make && make install
        if [ -d "tools" ]
        then
            cd "tools"
        fi
        python setup.py install --user &> /dev/null
        exit
        ;;
    c)
	    echo -e "${TC}------ create pre-compiled POLARIS ------${NC}"
        cd ${install_directory}
        # Check for necessary programms
        if ! command_exists cmake
        then
            echo -e "--- ${RED}Error:${NC} cmake not found. Please install cmake!"
            exit
        fi
        # Installing fits libraries
        install_fits_support
        install_polaris 'ON'
        exit
        ;;
    D)
        printf '%s\n' "Do you really want to delete your POLARIS installation [y/N]?"
        read really_delete
        case ${really_delete:=n} in
            [yY]*) 
                echo -e "${TC}------ delete POLARIS ------${NC}"
                export_str="export PATH=\"${install_directory}/bin:"'$PATH'"\""
                if grep -q "${export_str}" ${HOME}/.bashrc
                then
                    sed -i.bak "/${export_str//\//\\/}/d" ${HOME}/.bashrc
                fi
                export_str="export LD_LIBRARY_PATH=\"${install_directory}/CCfits/.libs/:${install_directory}/cfitsio:"'${LD_LIBRARY_PATH}'"\""
                if grep -q "${export_str}" ${HOME}/.bashrc
                then
                    sed -i.bak "/${export_str//\//\\/}/d" ${HOME}/.bashrc
                fi
                cd ${install_directory}/../
                rm -rv ${install_directory}
                pip uninstall PolarisTools
                exit
                ;;
            *) 
                exit
                ;; 
        esac
	    ;;
    \?)
        usage
        ;;
    esac
done

# ---------------------------------------------------------------------------------
# -------------------------- Installation of POLARIS ------------------------------
# ---------------------------------------------------------------------------------
echo -e "${TC}------ Installer for the radiative transfer code POLARIS ------${NC}"

# Go to install directory
cd ${install_directory}

# Make shell scripts executable (should not be necessary)
chmod +x *.sh

if [ ! -d "projects" ]
then
    mkdir "projects"
fi

# Ask for additional features
echo -e "${PC}--- Additional features ---${NC}"

printf '%s\n' "Do you want to enable PolarisTools [y/N]? (Python scripts collection)"
read install_tools
case ${install_tools:=n} in
    [yY]*) 
        polaris_tools=true
    ;;
    *) 
        polaris_tools=false
    ;; 
esac

# Install Polaris
echo -e "${PC}--- Install POLARIS ---${NC}"
if ${install_directory}/src/polaris > /dev/null
then 
    echo -e "POLARIS pre-compiled binary can be used [${GREEN}done${NC}]"
    cmake . > /dev/null
    make install > /dev/null
else
    install_fits_support
    install_polaris 'OFF'
fi

# install PolarisTools
if [ "${polaris_tools}" == true ]
then
    echo -e "${PC}--- Install PolarisTools ---${NC}"
    install_polaris_tools
    echo -ne "Setting up PolarisTools ... "\\r
    cd "tools/"
    python setup.py install --user > /dev/null 2>&1 \
        || { echo -e  "Setting up PolarisTools [${RED}Error${NC}]"; exit; }
fi  

echo -e  "${TC}-> Installation of POLARIS ${NC}[${GREEN}done${NC}]"

source ~/.bashrc

