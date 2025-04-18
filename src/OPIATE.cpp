/****************************************************************************/
/* Code:       POLARIS_v2                                                   */
/*                                                                          */
/* Author:     -Stefan Reissl                                               */
/*              reissl@uni-heidelberg.de                                    */
/*                                                                          */
/*              University of Heidelberg,                                   */
/*              Institute of Theoretical Astrophysics (ITA),                */
/*              Albert-Ueberle-Str. 2, 69120 Heidelberg, Germany            */
/*                                                                          */
/*                                                                          */
/* COPYRIGHT:                                                               */
/* The code is free of charge for any scientific purpose.                   */
/* This software is provided in the hope that it will be useful but         */
/* without any warranty of ability or fitness of a particular purpose.      */
/* We also reject any responsibility for incorrect results that may be      */
/* produced with this code.                                                 */
/* Any publication that makes use of the POLARIS software                   */
/* (completely or in part) must mention the name of the POLARIS code and    */
/* cite the paper Reissl et al. 2016                                        */
/* http://esoads.eso.org/abs/2016A%26A...593A..87R                          */
/*                                                                          */
/* History:   29.11.2016                                                    */
/****************************************************************************/

#include "OPIATE.hpp"

/*void OPIATE::formatLine(string &line)
{
    string::size_type pos = 0;

    if(line.find_first_of("#") != string::npos)
    {
        pos = line.find("#");
        line.erase(pos, line.length() - pos - 1);
    }

    if(line.find_first_of("!") != string::npos)
    {
        pos = line.find("!");
        line.erase(pos, line.length() - pos - 1);
    }

    if(line.size() < 3)
    {
        line = "";
        return;
    }

    while(line.find('\t') != string::npos)
    {
        pos = line.find('\t');
        line.replace(pos, 1, " ");
    }

    while(line.find(" \r\n") != string::npos)
    {
        pos = line.find(" \r\n");
        line.replace(pos, 3, " ");
    }

    while(line.find(" \r") != string::npos)
    {
        pos = line.find(" \r");
        line.replace(pos, 2, " ");
    }

    while(line.find(" \n") != string::npos)
    {
        pos = line.find(" \n");
        line.replace(pos, 2, " ");
    }

    while(line.find("\r\n") != string::npos)
    {
        pos = line.find("\r\n");
        line.replace(pos, 2, " ");
    }

    while(line.find("\r") != string::npos)
    {
        pos = line.find("\r");
        line.replace(pos, 1, " ");
    }

    while(line.find("\n") != string::npos)
    {
        pos = line.find("\n");
        line.replace(pos, 1, " ");
    }

    while(line.find("  ") != string::npos)
    {
        pos = line.find("  ");
        line.replace(pos, 2, " ");
    }

    if(line == " ")
        line = "";

    while(line.find(",") != string::npos)
    {
        pos = line.find(",");
        line.replace(pos, 1, ".");
    }

    if(line.size() > 0)
    {
        while(line.c_str()[line.size() - 1] == ' ')
        {
            pos = line.find_last_of(' ');
            line.erase(pos, 1);
        }
    }

    while(line.c_str()[0] == ' ')
    {
        line.erase(0, 1);
    }
}*/

void COpiateDataBase::printParameters(parameters & param, CGridBasic * grid)
{
    cout << CLR_LINE;
    cout << "OPIATE parameters                             " << endl;
    cout << SEP_LINE;

    if(database_counter==0)
    {
        cout << ERROR_LINE << "No OPIATE database available!               \n" ;
        return;
    }

    if(database_counter==1)
    {
        cout << "- Emission   database:\n   " << path_emi << endl;
        cout << "- Absorption database: none                 \n";
    }

    if(database_counter==2)
    {
        cout << "- Emission   database:\n   " << path_emi << endl;
        cout << "- Absorption database:\n   " << path_abs << endl;
    }

    cout << "- Velocity field     : ";

    if(grid->isVelocityFieldAvailable())
        cout << "taken from grid" << endl;
    else
        cout << "none" << endl;

    cout << "- Turbulent Velocity : ";
    if(param.getTurbulentVelocity() > 0)
        cout << param.getTurbulentVelocity() << " [m/s]" << endl;
    else if(grid->isTurbulentVelocityAvailable())
        cout << "taken from grid" << endl;
    else
        cout << "none" << endl;

    cout << "\n- Available species:    " << endl;
    for(uint i = 0; i<max_species; i++)
    {
        cout << "  - Species nr. " << i +1 << ", label: " <<  list_names[i];

        if(list_weight[i]>0)
            cout << ", weight: " << list_weight[i];

        double freq=list_freq[i];

        if(freq/1.0e9>1.0)
        {
            cout << ", freq.: " << freq/1.0e9 << " GHz" << endl;
        }
        else if(freq/1.0e6>1.0)
        {
            cout << ", freq.: " << freq/1.0e6 << " MHz" << endl;
        }
        else if(freq/1.0e3>1.0)
        {
            cout << ",freq.: " << freq/1.0e6 << " kHz" << endl;
        }
        else
        {
            cout << "   freq. : " << freq/1.0e6 << " Hz" << endl;
        }
    }

    cout << "\n- Selected species :    " << endl;

    for(uint i=0;i<param.getNrOfOPIATESpecies();i++)
    {
        cout << "  - Detector nr. " << i+1 << ", label: " << param.getOpiateSpec(i) << endl;
    }

    cout << SEP_LINE;
}


bool COpiateDataBase::readOpiateDataBase(parameters & param)
{
    string path_emi=param.getOpiatePathEmission();
    string path_abs=param.getOpiatePathAbsorption();

    if(path_emi.size()>0)
    {
        if(!readEmissivityData(path_emi))
            return false;

        has_emi_data=true;
    }
    else
    {
        cout << CLR_LINE;
        cout << ERROR_LINE << "A path to an OPIATE emissivity database is required!               \n" ;
        return false;
    }

    if(path_abs.size()>0)
    {
        if(!readAbsorptionData(path_abs))
            return false;

        has_abs_data=true;
    }

    return true;
}

bool COpiateDataBase::readDataBase(string filename)
{
//
    //cout << "Reading OPIATE fits data from:\n       " << filename <<  "               \n" << flush;
    //auto_ptr<FITS> pInfile(0);
    unique_ptr<FITS> pInfile;

    cout << CLR_LINE;
    cout << "-> Reading db nr.: " << database_counter+1 << " ...             \r" << flush;

    try
    {
        pInfile.reset(new FITS(filename.c_str(),Read,true));
    }
    catch(CCfits::FITS::CantOpen)
    {
        cout << CLR_LINE;
        cout << ERROR_LINE << "Cannot open OPIATE file:\n" << filename << "   \n" ;
        cout << "         Check path and file format!                   \n" ;
        return false;
    }

    PHDU& image = pInfile->pHDU();
    valarray<double>  contents;

    image.readAllKeys();
    image.read(contents);

    long max_row=image.axis(1);
    long max_col=image.axis(0);

    cout << contents.size() << endl << flush;
    uint max_species=uint((max_col-1)/2);

    if(database_counter==0)
    {
        max_ids=max_row;
        list_IDs=new double[max_ids];

        for(uint i=0;i<max_ids;i++)
            list_IDs[i]=contents[i];
    }
    else
    {
        //check with prev. data bases
    }

    database_counter++;

    cout << CLR_LINE;

    for(uint i=0; i<max_species;i++)
    {
        char str_tmp[32];
        char str_end[32];

        cout << "-> Reading db nr.: " << database_counter+1 << ", species: " << i+1 << "         \r" << flush;

        //copy for WINDOWS has to be adjusted here

        strcpy(str_tmp, "%03d");
        sprintf(str_end, str_tmp, i+1);

        string key_name="SNA_";
        string key_weight="SWE_";
        string key_freq="SFR_";

        key_name+=str_end;
        key_weight+=str_end;
        key_freq+=str_end;

        string s_name;
        double s_weight;
        double s_freq;

        try
        {
            image.readKey(key_name, s_name);
        }
        catch(CCfits::HDU::NoSuchKeyword)
        {
            cout << CLR_LINE;
            cout << ERROR_LINE << "Keyword \""<< key_name << "\" is required in file:\n       " << filename <<  "               \n" ;
            return false;
        }

        try
        {
            image.readKey(key_weight, s_weight);
        }
        catch(CCfits::HDU::NoSuchKeyword)
        {
            cout << CLR_LINE;
            cout << ERROR_LINE << "Keyword \""<< key_weight << "\" is required in file:\n       " << filename <<  "               \n" ;
            return false;
        }

        try
        {
            image.readKey(key_freq, s_freq);
        }
        catch(CCfits::HDU::NoSuchKeyword)
        {
            cout << CLR_LINE;
            cout << ERROR_LINE << "Keyword \""<< key_freq << "\" is required in file:\n       " << filename <<  "               \n" ;
            return false;
        }

        COpiateEntry * entry=new COpiateEntry(max_ids);

        entry->weight=s_weight;
        entry->freq=s_freq;
        entry->name=s_name;

        for(uint j=0;j<max_ids;j++)
        {
            double em=0;
            double ex=0;

            entry->setData(j,em,ex);
        }

        entries.push_back(entry);
    }



//    for(uint j=0;j<max_)
//
//    i * m_n + j


    //continue here ---
    cout << CLR_LINE;
    return true;
}

bool COpiateDataBase::readFitsData(string filename, Matrix2D & mat)
{
        //
        //cout << "Reading OPIATE fits data from:\n       " << filename <<  "               \n" << flush;
        // auto_ptr<FITS> pInfile(0);
        unique_ptr<FITS> pInfile;

        cout << CLR_LINE;
        cout << "-> Reading fits data...              \r" << flush;

        try
        {
            pInfile.reset(new FITS(filename.c_str(),Read,true));
        }
        catch(CCfits::FITS::CantOpen)
        {
            cout << CLR_LINE;
            cout << ERROR_LINE << "Cannot open OPIATE file:\n" << filename << "   \n" ;
            cout << "         Check path and file format!                   \n" ;
            return false;
        }

        PHDU& image = pInfile->pHDU();
        valarray<double>  contents;

        image.readAllKeys();
        image.read(contents);

        long max_row=image.axis(1);
        long max_col=image.axis(0);

        if(database_counter==0)
        {
            max_species=max_col-1;
            max_ids=max_row;

            list_freq=new double[max_species];
            list_weight=new double[max_species];
            list_IDs=new double[max_row];
        }

        if(database_counter==1)
        {
            if(max_species!=max_col-1)
            {
                cout << CLR_LINE;
                cout << ERROR_LINE << "Ammount of species in differen OPIATE files do not match!            \n" ;
                cout << "         Identical order and values are required in all files!           \n" ;
                return false;
            }

            if(max_ids!=max_row)
            {
                cout << CLR_LINE;
                cout << ERROR_LINE << "Ammount of unique OPIATE files do not match!            \n" ;
                cout << "         Identical order and values are required in all files!           \n" ;
                return false;
            }
        }

        if(database_counter>1)
        {
            cout << CLR_LINE;
            cout << ERROR_LINE << "Only two databases are currently supported!              \n" ;
            return false;
        }

        cout << CLR_LINE;
        cout << "-> Scanning for keys ...              \r" << flush;

        for(uint i=0; i<max_species;i++)
        {
            char str_tmp[1024];
            char str_end[1024];

            //copy for WINDOWS has to be adjusted here

            strcpy(str_tmp, "%03d");
            sprintf(str_end, str_tmp, i+1);

            string key_name="SNA_";
            string key_weight="SWE_";
            string key_freq="SFR_";

            key_name+=str_end;
            key_weight+=str_end;
            key_freq+=str_end;

            string s_name;
            double s_weight;
            double s_frequ;

            try
            {
                image.readKey(key_name, s_name);
            }
            catch(CCfits::HDU::NoSuchKeyword)
            {
                cout << CLR_LINE;
                cout << ERROR_LINE << "Keyword \""<< key_name << "\" is required in file:\n       " << filename <<  "               \n" ;
                return false;
            }

            try
            {
                image.readKey(key_weight, s_weight);
            }
            catch(CCfits::HDU::NoSuchKeyword)
            {
                cout << CLR_LINE;
                cout << ERROR_LINE << "Keyword \""<< key_weight << "\" is required in file:\n       " << filename <<  "               \n" ;
                return false;
            }

            try
            {
                image.readKey(key_freq, s_frequ);
            }
            catch(CCfits::HDU::NoSuchKeyword)
            {
                cout << CLR_LINE;
                cout << ERROR_LINE << "Keyword \""<< key_freq << "\" is required in file:\n       " << filename <<  "               \n" ;
                return false;
            }

            //cout << "Data:" << key_name << ":\t\"" << s_name << "\"\n" << flush;
            //cout << "Data:" << key_weight << ":\t\"" << s_weight << "\"\n" << flush;
            //cout << "Data:" << key_freq << ":\t\"" << s_frequ << "\"\n" << flush;

            if(database_counter==0)
            {
                list_names.push_back(s_name);
                list_freq[i]=s_frequ;
                list_weight[i]=s_weight;
            }
            else
            {
                if(list_freq[i]!=s_frequ)
                {
                    cout << CLR_LINE;
                    cout << ERROR_LINE << "Frequencies of different OPIATE files do not match!            \n" ;
                    cout << "         Identical order and values are required in all files!           \n" ;
                    return false;
                }

                if(list_names[i].compare(s_name)!=0)
                {
                    cout << CLR_LINE;
                    cout << ERROR_LINE << "Names of different OPIATE files do not match!            \n" ;
                    cout << "         Identical order and values are required in all files!           \n" ;
                    return false;
                }

                if(list_weight[i]!=s_weight)
                {
                    cout << CLR_LINE;
                    cout << ERROR_LINE << "Weigths of different OPIATE files do not match!            \n" ;
                    cout << "         Identical order and values are required in all files!           \n" ;
                    return false;
                }


            }

            //tmp_mat.set(i,contents[i]);

        }

        Matrix2D tmp_mat;

        //cout << "row: " << max_row << "\n";
        //cout << "col: " << max_col << "\n";

        mat.resize(max_row,max_col-1);
        tmp_mat.resize(max_row,max_col);
        //cout << contents.size() << "\n";

        cout << CLR_LINE;
        cout << "Creating matrix ...                \r" << flush;
        for (long i = 0; i < max_col*max_row; i++)
        {
            tmp_mat.set(i,contents[i]);
        }

        //cout << "\n";

        //tmp_mat.printMatrix();

        cout << CLR_LINE;
        cout << "-> Creating ID table ...              \r" << flush;
        for (long i = 0; i < max_row; i++)
        {
            if(database_counter==0)
            {
                list_IDs[i]=tmp_mat(i,0);
            }
            else
            {
                if(list_IDs[i]!=tmp_mat(i,0))
                {
                    cout << ERROR_LINE << "Unique OPIATE IDs of different OPIATE files do not match!            \n" ;
                    return false;
                }
            }

            //cout << i << "\t" <<list_IDs[i] << "\n";

            if(i>0)
            {
                if(list_IDs[i-1]>list_IDs[i])
                {
                    cout << CLR_LINE;
                    cout << ERROR_LINE << "OPIATE IDs are not in ascending order!               \n" ;
                    cout << "         Check first column of your input fits file!                   \n" ;
                    return false;
                }

                if(list_IDs[i-1]==list_IDs[i])
                {
                    cout << CLR_LINE;
                    cout << ERROR_LINE << "Identical OPIATE IDs detected!               \n" ;
                    cout << "         Check first column of your input fits file!                   \n" ;
                    return false;
                }
            }
        }

        cout << "\nMatrix:\n";

        for (long i = 0; i < max_row; i++)
        {
            for (long j = 1; j < max_col; j++)
            {
                mat(i,j-1) = tmp_mat(i,j);
            }
        }

        //cout << "\n";

        //mat.printMatrix();


    database_counter++;
    return true;
}

double COpiateDataBase::getGaussA(double temp_gas, double v_turb)
{
    double v_th = sqrt(2.0 * con_kB * temp_gas / (list_weight[current_index] * 1e-3 / con_Na));
    double gauss_a = 1.0 / sqrt(pow(v_th, 2) + pow(v_turb, 2));
    return gauss_a;
}

void COpiateDataBase::calcLineBroadening(CGridBasic * grid)
{
    long max_cells = grid->getMaxDataCells();
    
    cout << CLR_LINE;
    cout << "-> Calculating line broadening for each grid cell ...     \r" << flush;
    
    #pragma omp parallel for schedule(dynamic)
    for(long i_cell = 0; i_cell < long(max_cells); i_cell++)
    {
        cell_basic * cell = grid->getCellFromIndex(i_cell);

        // Get necessary quantities from the current cell
        double temp_gas = grid->getGasTemperature(*cell);
        double turbulent_velocity = grid->getTurbulentVelocity(cell);

        // Set gauss_a for each transition only once
        grid->setGaussA(cell, getGaussA(temp_gas, turbulent_velocity));
    }
    
    cout << CLR_LINE;
}

double COpiateDataBase::getProjCellVelocityInterp(const Vector3D & tmp_pos,
                                    const Vector3D & dir_map_xyz,
                                    const VelFieldInterp & vel_field_interp)
{
    double cell_velocity = 0;

    // Get the velocity in the photon direction of the current position
    if(vel_field_interp.vel_field.size() > 0 && !vel_field_interp.zero_vel_field)
    {
        // Get velocity from grid cell with interpolation
        Vector3D rel_pos = tmp_pos - vel_field_interp.start_pos;
        cell_velocity = vel_field_interp.vel_field.getValue(rel_pos.length());
    }
    return cell_velocity;
}

uint COpiateDataBase::getIndex(uint op_id) const
{
    uint N = max_ids;
    uint min = 0, max = N - 1;

    if(op_id < list_IDs[0] || op_id > list_IDs[max])
        return MAX_UINT;

    while(max - min > 1)
    {
        uint i = min + (max - min) / 2;
        if(list_IDs[i] >= op_id)
            max = i;
        else
            min = i;
    }

    if(list_IDs[min] == op_id)
        return min;

    uint lower = min - 1;
    uint upper = min + 1;

    if(lower == MAX_UINT)
        lower = 0;

    if(upper >= max_ids)
        upper = max_ids - 1;

    if(list_IDs[lower] == op_id)
        return lower;

    if(list_IDs[upper] == op_id)
        return upper;

    return min;
}

void COpiateDataBase::getMatrices(CGridBasic * grid,
                                const photon_package * pp,
                                uint i_spec,
                                uint i_trans,
                                double velocity,
                                const LineBroadening & line_broadening,
                                const MagFieldInfo & mag_field_info,
                                StokesVector * line_emissivity,
                                Matrix2D * line_absorption_matrix) const
{
    double emission = 0.0;
    double absorption = 0.0;
    
    uint op_id=grid->getOpiateID(pp);
    uint index = getIndex(op_id);
    
    // Reset absorption matrix and emissivity
    line_absorption_matrix->resize(4, 4);
    *line_emissivity = 0;
            
    if(index!=MAX_UINT)
    {
        if(has_abs_data==true)
            absorption = mat_absorption(index,current_index)*line_broadening.gauss_a;

        if(has_emi_data==true)
            emission = mat_emissivity(index,current_index)*line_broadening.gauss_a;
    }

    // Calculate the line matrix from rotation polarization matrix and line shape
    // getGaussLineMatrix(grid, pp, velocity, line_absorption_matrix);
    double line_amplitude = exp(-(pow(velocity, 2) * pow(line_broadening.gauss_a, 2))) / PIsq;
    
    // Only diagonal without polarization rotation matrix elements
    for(uint i = 0; i < 4; i++)
    {
        line_absorption_matrix->setValue(i, i, line_amplitude);
    }

    // Calculate the Emissivity of the gas particles in the current cell
    *line_emissivity = *line_absorption_matrix * StokesVector(emission, 0, 0, 0);

    // Calculate the line matrix from rotation polarization matrix and line shape
    *line_absorption_matrix *= absorption;
}

uint COpiateDataBase::getMaxSpecies()
{
    return max_species;
}

double COpiateDataBase::getFrequency(uint pos)
{
    return list_freq[pos];
}

bool COpiateDataBase::findIndexByName(string name)
{
    for(uint i = 0; i < list_names.size(); i++)
    {
        if(name.compare(list_names[i]) == 0)
        {
            current_index = i;
            return true;
        }
    }

    cout << CLR_LINE;
    cout << ERROR_LINE << "A species by the name of \"" << name
            << "\" is not listed in the loaded OPIATE databases!              \n";
    return false;
}

string COpiateDataBase::getCurrentName()
{
    return list_names[current_index];
}

bool COpiateDataBase::readEmissivityData(string filename)
{
    path_emi=filename;
    return readFitsData(filename, mat_emissivity);
}

bool COpiateDataBase::readAbsorptionData(string filename)
{
    path_abs=filename;
    return readFitsData(filename, mat_absorption);
}

double COpiateDataBase::getCurrentFrequency()
{
    return list_freq[current_index];
}

double COpiateDataBase::getMolecularWeight()
{
    return list_weight[current_index];
}

bool COpiateDataBase::initVelChannels(uint nr_of_channels, double max_vel)
{
    if(velocity_channel != 0)
    {
        delete[] velocity_channel;
        velocity_channel = 0;
    }

    nr_of_velocity_channels = nr_of_channels;
    max_velocity = max_vel;

    velocity_channel = new double[nr_of_velocity_channels];

    if(nr_of_velocity_channels > 1)
    {
        for(uint i = 0; i < nr_of_velocity_channels; i++)
        {
            velocity_channel[i] =
                2 * (float)i / ((float)nr_of_velocity_channels - 1) * max_velocity - max_velocity;
        };
    }
    else if(nr_of_velocity_channels == 1)
    {
        velocity_channel[0] = 0;
    }
    else
    {
        cout << CLR_LINE;
        cout << ERROR_LINE << "Number of velocity channels is not larger than zero!                 \n";
        return false;
    }

    return true;
}

double COpiateDataBase::getVelocityChannel(uint vch)
{
    return velocity_channel[vch];
}

uint COpiateDataBase::getNrOfVelChannels()
{
    return nr_of_velocity_channels;
}

double COpiateDataBase::getMaxVelocity()
{
    return max_velocity;
}
