/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#include <stdlib.h>
#include "GridCylindrical.hpp"

bool CGridCylindrical::loadGridFromBinaryFile(parameters & param, uint _data_len)
{
    ushort tmpID, tmpOffset;
    string filename = param.getPathGrid();
    bool b_center = false;

    uint r_counter = 0;
    uint ph_counter = 0;
    uint z_counter = 0;

    line_counter = 0;
    char_counter = 0;

    ifstream bin_reader(filename.c_str(), ios::in | ios::binary);

    if(bin_reader.fail())
    {
        cout << ERROR_LINE << "Cannot load binary cylindrical grid file:" << endl;
        cout << filename << "\n" << endl;
        return false;
    }

    resetGridValues();

    turbulent_velocity = param.getTurbulentVelocity();

    min_len = 1e300;

    max_cells = 0;

    line_counter = 1;
    char_counter = 0;
    float last_percentage = 0;

    bin_reader.read((char *)&tmpID, 2);
    bin_reader.read((char *)&tmpOffset, 2);

    dataID = tmpID;
    data_offset = (uint)tmpOffset;
    data_len = _data_len + data_offset;

    if(dataID == GRID_ID_CYL)
    {
        data_ids.resize(data_offset);

        for(uint i = 0; i < data_offset; i++)
        {
            ushort tmp_ids = 0;
            bin_reader.read((char *)&tmp_ids, 2);
            data_ids[i] = tmp_ids;
        }

        if(!setDataPositionsVariable())
            return false;
    }
    else
    {
        cout << ERROR_LINE << "A cylindrical grid requires an ID of \"" << GRID_ID_CYL << "\"!" << endl;
        return false;
    }

    uint tmp_data_offset = validateDataPositions(param);
    if(tmp_data_offset == MAX_UINT)
        return false;

    bin_reader.read((char *)&Rmin, 8);
    bin_reader.read((char *)&Rmax, 8);
    bin_reader.read((char *)&Zmax, 8);
    bin_reader.read((char *)&N_r, 2);

    // Get global number of phi cells
    uint N_ph_tmp = 0;
    bin_reader.read((char *)&N_ph_tmp, 2);

    // Set number of phi cells in each ring
    N_ph = new uint[N_r];
    for(uint i_r = 0; i_r < N_r; i_r++)
        N_ph[i_r] = N_ph_tmp;

    bin_reader.read((char *)&N_z, 2);
    bin_reader.read((char *)&log_factorR, 8);
    bin_reader.read((char *)&log_factorPh, 8);
    bin_reader.read((char *)&log_factorZ, 8);

    // Convert limits with conversion factor
    Rmin *= conv_length_in_SI;
    Rmax *= conv_length_in_SI;
    Zmax *= conv_length_in_SI;

    // --------------------------------------
    // ---------- Radial-direction ----------
    // --------------------------------------

    // Init radial cell border
    listR = new double[N_r + 1];
    if(log_factorR == 0)
    {
        // Allow user defined radius list, if log_factorR is zero

        // The global borders are already in the grid
        listR[0] = Rmin;
        listR[N_r] = Rmax;

        // Set the cell borders
        for(uint i_r = 1; i_r < N_r; i_r++)
        {
            // Read radial cell border position
            bin_reader.read((char *)&listR[i_r], 8);

            // Update radial position with conversion factors
            listR[i_r] *= conv_length_in_SI;
        }
    }
    else if(log_factorR == 1.0)
    {
        // Sinus shaped list, which emphasizes the middle rings
        CMathFunctions::SinList(Rmin, Rmax, listR, N_r + 1, log_factorR);
    }
    else if(log_factorR > 1.0)
    {
        // Exponentially increasing width of the cells in radial direction
        CMathFunctions::ExpList(Rmin, Rmax, listR, N_r + 1, log_factorR);
    }
    else
    {
        // Linear width of the cells in radial direction
        CMathFunctions::LinearList(Rmin, Rmax, listR, N_r + 1);
    }

    // -----------------------------------
    // ---------- Phi-direction ----------
    // -----------------------------------

    // Init phi cell border
    listPh = new double *[N_r];
    if(log_factorPh == 0)
    {
        // Allow user defined phi list, if log_factorPh is zero
        for(uint i_r = 0; i_r < N_r; i_r++)
        {
            // Init 2D cell borders in phi-direction
            listPh[i_r] = new double[N_ph[0] + 1];

            // The global borders are already in the grid
            listPh[i_r][0] = 0;
            listPh[i_r][N_ph[0]] = PIx2;
        }

        for(uint i_ph = 1; i_ph < N_ph[0]; i_ph++)
        {
            // Read cell border in z-direction
            double ph = 0;
            bin_reader.read((char *)&ph, 8);

            // Set the cell borders
            for(uint i_r = 0; i_r < N_r; i_r++)
                listPh[i_r][i_ph] = ph;
        }
    }
    else if(log_factorPh == -1.0)
    {
        // Allow user defined z list based on a variable dz, if log_factorZ is -1
        for(uint i_r = 0; i_r < N_r; i_r++)
        {
            // Read number of phi cells in the current ring
            bin_reader.read((char *)&N_ph[i_r], 2);

            // Init 2D cell borders in z-direction
            listPh[i_r] = new double[N_ph[i_r] + 1];
            CMathFunctions::LinearList(0, PIx2, listPh[i_r], N_ph[i_r] + 1);
        }
    }
    else
    {
        // Linear width of the cells in z-direction
        for(uint i_r = 0; i_r < N_r; i_r++)
        {
            listPh[i_r] = new double[N_ph[i_r] + 1];
            CMathFunctions::LinearList(0, PIx2, listPh[i_r], N_ph[i_r] + 1);
        }
    }

    // ---------------------------------
    // ---------- Z-direction ----------
    // ---------------------------------

    // Init z cell border
    listZ = new double *[N_r];
    if(log_factorZ == 0)
    {
        // Allow user defined z list, if log_factorZ is zero
        for(uint i_r = 0; i_r < N_r; i_r++)
        {
            // Init 2D cell borders in z-direction
            listZ[i_r] = new double[N_z + 1];

            // The global borders are already in the grid
            listZ[i_r][0] = -Zmax;
            listZ[i_r][N_z] = Zmax;
        }

        for(uint i_z = 1; i_z < N_z; i_z++)
        {
            // Read cell border in z-direction
            double z = 0;
            bin_reader.read((char *)&z, 8);

            // Set the cell borders
            for(uint i_r = 0; i_r < N_r; i_r++)
                listZ[i_r][i_z] = z * conv_length_in_SI;
        }
    }
    else if(log_factorZ == -1.0)
    {
        // Add an empty cell at both ends
        N_z += 2;

        // Allow user defined z list based on a variable dz, if log_factorZ is -1
        for(uint i_r = 0; i_r < N_r; i_r++)
        {
            // Init 2D cell borders in z-direction
            listZ[i_r] = new double[N_z + 1];

            // The global borders are already in the grid
            listZ[i_r][0] = -Zmax;
            listZ[i_r][N_z] = Zmax;

            // Read distance between two borders in z-direction for the current radius
            // cell
            double dz = 0;
            bin_reader.read((char *)&dz, 8);

            // Calculate the maximum height up to which information is in the grid
            double local_zmax = (N_z - 2) / 2. * dz * conv_length_in_SI;

            // Set the second and second last border from dz
            listZ[i_r][1] = -local_zmax;
            listZ[i_r][N_z - 1] = local_zmax;

            // Set the cell borders
            for(uint i_z = 2; i_z <= N_z - 2; i_z++)
                listZ[i_r][i_z] = -local_zmax + (i_z - 1) * dz;
        }
    }
    else if(log_factorZ == 1.0)
    {
        // Sinus shaped list, which emphasizes the midplane
        for(uint i_r = 0; i_r < N_r; i_r++)
        {
            listZ[i_r] = new double[N_z + 1];
            CMathFunctions::SinList(-Zmax, Zmax, listZ[i_r], N_z + 1, log_factorZ);
        }
    }
    else if(log_factorZ > 1.0)
    {
        // Exponentially increasing width of the cells in z-direction (symmetrically)
        for(uint i_r = 0; i_r < N_r; i_r++)
        {
            listZ[i_r] = new double[N_z + 1];
            CMathFunctions::ExpListSym(-Zmax, Zmax, listZ[i_r], N_z + 1, log_factorZ);
        }
    }
    else
    {
        // Linear width of the cells in z-direction
        for(uint i_r = 0; i_r < N_r; i_r++)
        {
            listZ[i_r] = new double[N_z + 1];
            CMathFunctions::LinearList(-Zmax, Zmax, listZ[i_r], N_z + 1);
        }
    }

    // -----------------------------------------
    // ---------- Check of the limits ----------
    // -----------------------------------------
    if(Rmin <= 0)
    {
        cout << ERROR_LINE << "Inner radius (Rmin = " << Rmin << " must be larger than zero!" << endl;
        return false;
    }

    if(Rmax <= 0)
    {
        cout << ERROR_LINE << "Outer radius (Rmax = " << Rmax << " must be larger than zero!" << endl;
        return false;
    }

    if(Rmax <= Rmin)
    {
        cout << ERROR_LINE << "Outer radius (Rmax = " << Rmax
             << ") must be larger than inner radius (Rmin = " << Rmin << ")!" << endl;
        return false;
    }

    if(Zmax <= 0)
    {
        cout << ERROR_LINE << "Maximum vertical extent must be larger than zero (Zmax = " << Zmax << ")!" << endl;
        return false;
    }

    if(N_r < 3)
    {
        cout << ERROR_LINE << "Number of cells in R direction has to be larger than 2!" << endl;
        return false;
    }

    // Init grid cells
    grid_cells = new cell_cyl ***[N_r];

    // Init center cells
    center_cells = new cell_cyl *[N_z];

    for(uint i_r = 0; i_r < N_r; i_r++)
    {
        grid_cells[i_r] = new cell_cyl **[N_ph[i_r]];

        // cout << "Allocating memory for cylindrical grid cells: " << float(100.0 * double(i_r) / double(N_r))
        //      << "      \r" << flush;

        for(uint i_ph = 0; i_ph < N_ph[i_r]; i_ph++)
        {
            grid_cells[i_r][i_ph] = new cell_cyl *[N_z];

            for(uint i_z = 0; i_z < N_z; i_z++)
            {
                grid_cells[i_r][i_ph][i_z] = 0;
                center_cells[i_z] = 0;
            }
        }
    }

    // Clear user output
    cout << CLR_LINE;

    max_cells = N_z;
    for(uint i_r = 0; i_r < N_r; i_r++)
        max_cells += N_ph[i_r] * N_z;

    line_counter = -1;

    max_len = max(2 * Rmax, 2 * Zmax);
    min_len = listR[1] - listR[0];

    total_volume = PI * Rmax * Rmax * (2 * Zmax);

    while(!bin_reader.eof())
    {
        line_counter++;

        if(r_counter < N_r)
        {
            double dr = listR[r_counter + 1] - listR[r_counter];

            if(dr == 0)
            {
                cout << ERROR_LINE << "No step size in r-direction of cylindrical grid!" << endl;
                return false;
            }

            if(dr < min_len)
                min_len = dr;

            if(ph_counter < N_ph[r_counter])
            {
                double dph = listPh[r_counter][ph_counter + 1] - listPh[r_counter][ph_counter];

                if(dph == 0)
                {
                    cout << ERROR_LINE << "No step size in phi-direction of cylindrical grid!" << endl;
                    return false;
                }

                double d = 2 * Rmin * dph / PIx2;

                if(d < min_len)
                    min_len = d;
            }
        }

        if(r_counter == 0)
        {
            if(z_counter < N_z)
            {
                double dz = listZ[0][z_counter + 1] - listZ[0][z_counter];

                if(dz == 0)
                {
                    cout << ERROR_LINE << "No step size in z-direction of cylindrical grid!" << endl;
                    cout << "Please look in the manual" << endl;
                    return false;
                }

                if(dz < min_len)
                    min_len = dz;
            }
        }

        if(z_counter == N_z && b_center == false)
        {
            ph_counter++;
            z_counter = 0;
        }

        if(r_counter < N_r && ph_counter == N_ph[r_counter] && b_center == false)
        {
            r_counter++;
            ph_counter = 0;
        }

        if(r_counter == N_r && b_center == false)
        {
            b_center = true;
            z_counter = 0;
        }

        if(z_counter == N_z && b_center == true)
            break;

        // Calculate percentage of total progress per source
        float percentage = 100.0 * double(line_counter) / double(max_cells);

        // Show only new percentage number if it changed
        if((percentage - last_percentage) > PERCENTAGE_STEP)
        {
            char_counter++;
            cout << "-> Loading cylindrical grid file: " << percentage << " [%]      \r" << flush;
            last_percentage = percentage;
        }

        cell_cyl * tmp_cell = 0;

        if(b_center == true)
        {
            center_cells[z_counter] = new cell_cyl;
            center_cells[z_counter]->setRID(MAX_UINT);
            center_cells[z_counter]->setPhID(MAX_UINT);
            center_cells[z_counter]->setZID(z_counter);
            center_cells[z_counter]->resize(data_len + tmp_data_offset);
            tmp_cell = center_cells[z_counter];
            tmp_cell->setID(line_counter);
        }
        else
        {
            grid_cells[r_counter][ph_counter][z_counter] = new cell_cyl;
            grid_cells[r_counter][ph_counter][z_counter]->setRID(r_counter);
            grid_cells[r_counter][ph_counter][z_counter]->setPhID(ph_counter);
            grid_cells[r_counter][ph_counter][z_counter]->setZID(z_counter);
            grid_cells[r_counter][ph_counter][z_counter]->resize(data_len + tmp_data_offset);
            tmp_cell = grid_cells[r_counter][ph_counter][z_counter];
        }

        if(log_factorZ != -1 || (z_counter < N_z - 1 && z_counter > 0))
            for(uint i = 0; i < data_offset; i++)
            {
                double tmp_data1 = 0;
                bin_reader.read((char *)&tmp_data1, 8);
                tmp_cell->setData(i, tmp_data1);
            }

        updateVelocity(tmp_cell, param);

        if(uint(tmp_cell->getData(data_pos_id)) < 0 ||
           uint(tmp_cell->getData(data_pos_id)) > param.getMaxDustComponentChoice())
        {
            cout << ERROR_LINE << "Dust ID in grid exceeds maximum number of dust choices "
                    "available! "
                 << endl;
            return false;
        }

        tmp_cell->setID(line_counter);

        updateDataRange(tmp_cell);

        double tmp_vol = getVolume(*tmp_cell);
        total_gas_mass += getGasMassDensity(*tmp_cell) * tmp_vol;
        cell_volume += tmp_vol;
        z_counter++;
    }

    bin_reader.close();

    if(max_cells != uint(line_counter))
    {
        cout << ERROR_LINE << "Number of read in cells do not match the maximal number of "
                "expected cells!"
             << endl;
        cout << "       Expected " << max_cells << " cells, but found " << uint(line_counter)
             << " cells in grid!" << endl;
        return false;
    }

    if(min_len < 0)
    {
        cout << ERROR_LINE << "Minimum length is smaller than zero!" << endl;
        return false;
    }

    data_len += tmp_data_offset;
    data_offset += tmp_data_offset;

    cout << CLR_LINE;
    return true;
}

bool CGridCylindrical::writePlotFiles(string path, parameters & param)
{
    nrOfPlotPoints = param.getNrOfPlotPoints();
    nrOfPlotVectors = param.getNrOfPlotVectors();
    maxPlotLines = param.getMaxPlotLines();

    if(nrOfPlotPoints + nrOfPlotVectors == 0)
        return true;

    if(max_cells == 0)
    {
        cout << ERROR_LINE << "Cannot write cylindrical plot file to:" << endl;
        cout << path;
        cout << "Not enough tree cells available! " << endl;
        return false;
    }

    plt_gas_dens = (size_gd_list > 0);  // 1
    plt_mol_dens = (nrOfDensRatios>0);
    plt_dust_dens = false;                       // param.getPlot(plIDnd) && (!data_pos_dd_list.empty()); // 2
    plt_gas_temp = (data_pos_tg != MAX_UINT);    // 3
    plt_dust_temp = (!data_pos_dt_list.empty()); // 4
    plt_rat = (!data_pos_aalg_list.empty());     // 5
    plt_delta = (data_pos_tg != MAX_UINT) && (data_pos_mx != MAX_UINT); // && (data_pos_td != MAX_UINT); // 6
    plt_larm = (data_pos_tg != MAX_UINT) && (data_pos_mx != MAX_UINT);  // && (data_pos_td != MAX_UINT); // 7
    plt_mach = (data_pos_vx != MAX_UINT) && (data_pos_tg != MAX_UINT);  // 8

    plt_mag = (data_pos_mx != MAX_UINT); // 0
    plt_vel = (data_pos_vx != MAX_UINT); // 1

    if(nrOfPlotPoints <= 1)
    {
        nrOfPlotPoints = max_cells / 10;

        plt_gas_dens = false;
        plt_dust_dens = false;
        plt_gas_temp = false;
        plt_dust_temp = false;
        plt_rat = false;
        plt_delta = false;
        plt_larm = false;
        plt_mach = false;
    }
    else
        nrOfPlotPoints = max_cells / nrOfPlotPoints;

    if(nrOfPlotVectors <= 1)
    {
        nrOfPlotVectors = max_cells / 10;
        plt_mag = false;
        plt_vel = false;
    }
    else
        nrOfPlotVectors = max_cells / nrOfPlotVectors;

    if(nrOfPlotPoints == 0)
        nrOfPlotPoints = 1;

    if(nrOfPlotVectors == 0)
        nrOfPlotVectors = 1;

    stringstream point_header, vec_header, basic_grid_l0, basic_grid_l1;

    string grid_filename = path + "grid_geometry.plt";
    string dens_gas_filename = path + "grid_gas_density.plt";
    string dens_dust_filename = path + "grid_dust_density.plt";
    string temp_gas_filename = path + "grid_gas_temp.plt";
    string temp_dust_filename = path + "grid_dust_temp.plt";
    string rat_filename = path + "grid_RAT.plt";
    string delta_filename = path + "grid_data.dat";
    string larm_filename = path + "grid_mag.plt";
    string mach_filename = path + "grid_vel.plt";
    string mag_filename = path + "grid_mag.plt";
    string vel_filename = path + "grid_vel.plt";

    ofstream point_fields[9];
    ofstream vec_fields[2];

    point_fields[0].open(grid_filename.c_str());

    if(point_fields[0].fail())
    {
        cout << ERROR_LINE << "Cannot write to:\n " << grid_filename << endl;
        return false;
    }

    if(plt_gas_dens)
    {
        point_fields[1].open(dens_gas_filename.c_str());

        if(point_fields[1].fail())
        {
            cout << ERROR_LINE << "Cannot write to:\n " << dens_gas_filename << endl;
            return false;
        }
    }

    if(plt_dust_dens)
    {
        point_fields[2].open(dens_dust_filename.c_str());

        if(point_fields[2].fail())
        {
            cout << ERROR_LINE << "Cannot write to:\n " << dens_dust_filename << endl;
            return false;
        }
    }

    if(plt_gas_temp)
    {
        point_fields[3].open(temp_gas_filename.c_str());

        if(point_fields[3].fail())
        {
            cout << ERROR_LINE << "Cannot write to:\n " << temp_gas_filename << endl;
            return false;
        }
    }

    if(plt_dust_temp)
    {
        point_fields[4].open(temp_dust_filename.c_str());

        if(point_fields[4].fail())
        {
            cout << ERROR_LINE << "Cannot write to:\n " << temp_dust_filename << endl;
            return false;
        }
    }

    if(plt_rat)
    {
        point_fields[5].open(rat_filename.c_str());

        if(point_fields[5].fail())
        {
            cout << ERROR_LINE << "Cannot write to:\n " << rat_filename << endl;
            return false;
        }
    }

    if(plt_delta)
    {
        point_fields[6].open(delta_filename.c_str());

        if(point_fields[6].fail())
        {
            cout << ERROR_LINE << "Cannot write to:\n " << delta_filename << endl;
            return false;
        }
    }

    if(plt_larm)
    {
        point_fields[7].open(larm_filename.c_str());

        if(point_fields[7].fail())
        {
            cout << ERROR_LINE << "Cannot write to:\n " << larm_filename << endl;
            return false;
        }
    }

    if(plt_mach)
    {
        point_fields[8].open(mach_filename.c_str());

        if(point_fields[8].fail())
        {
            cout << ERROR_LINE << "Cannot write to:\n " << mach_filename << endl;
            return false;
        }
    }

    if(plt_mag)
    {
        vec_fields[0].open(mag_filename.c_str());

        if(vec_fields[0].fail())
        {
            cout << ERROR_LINE << "Cannot write to:\n " << mag_filename << endl;
            return false;
        }
    }

    if(plt_vel)
    {
        vec_fields[1].open(vel_filename.c_str());

        if(vec_fields[1].fail())
        {
            cout << ERROR_LINE << "Cannot write to:\n " << vel_filename << endl;
            return false;
        }
    }

    line_counter = 0;
    char_counter = 0;

    point_header.str("");
    point_header << "reset" << endl;
    point_header << "#set terminal postscript" << endl;
    point_header << "#set output \'filename.plt\'" << endl;
    point_header << "set ticslevel 0" << endl;
    point_header << "set size ratio -1" << endl;
    point_header << "set view 45,45" << endl;

    point_header << "set xlabel \'x[m]\'" << endl;
    point_header << "set ylabel \'y[m]\'" << endl;
    point_header << "set zlabel \'z[m]\'" << endl;

    point_header << "set xrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;
    point_header << "set yrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;
    point_header << "set zrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;

    point_header << "set style arrow 1 nohead lt 2 lc rgb 0x0000ff lw 0.5" << endl;
    point_header << "set style arrow 2 nohead lt 2 lc rgb 0x0054dd lw 0.5" << endl;

    point_header << "set style line 1 pt 1 ps variable lt palette" << endl;

    point_header << "set grid" << endl;
    point_header << "set nokey" << endl;

    // 0 spherical grid
    point_fields[0] << point_header.str();
    point_fields[0] << "set title \'3D cylindrical grid geometry\' font \'Arial,12\'" << endl;
    point_fields[0] << "set style arrow 3 nohead ls 1 lw 0.5 lc rgb 0x550066" << endl;
    point_fields[0] << "splot '-' with vectors as 3,'-' with vectors as 2,'-' with vectors as 1" << endl;

    // 1 gas density
    point_fields[1] << point_header.str();
    point_fields[1] << "set title \'3D gas number density distribution (min: " << min_gas_dens
                    << "[m^-3]; max: " << max_gas_dens << "[m^-3])\' font \'Arial,12\'" << endl;
    point_fields[1] << "set cblabel \'gas density[m^-3]\'" << endl;
    point_fields[1] << "set palette defined (0 0.5 0 0, 1 0 0 1, 2 0 1 1, 3 1 1 0)" << endl;

    if(min_gas_dens == 0 && max_gas_dens == 0)
    {
        min_gas_dens = 0.1;
        max_gas_dens = 1;
    }

    if(min_gas_dens == 0)
    {
        min_gas_dens = 0.001 * max_gas_dens;
    }

    if(min_gas_dens / max_gas_dens > 0.9)
        min_gas_dens = 0.9 * max_gas_dens;

    point_fields[1] << "set cbrange[" << log10(min_gas_dens) << ":" << log10(max_gas_dens) << "]" << endl;
    point_fields[1] << "set format cb \'%.02g\'" << endl;

    point_fields[1] << "splot  '-' w p ls 1,'-' with vectors as 2,'-' with vectors as 1" << endl;

    // 2 dust density
    point_fields[2] << point_header.str();
    point_fields[2] << "set title \'3D gas number density distribution (min: " << min_dust_dens
                    << "[m^-3]; max: " << max_dust_dens << "[m^-3])\' font \'Arial,12\'" << endl;
    point_fields[2] << "set cblabel \'gas density[m^-3]\'" << endl;
    point_fields[2] << "set palette defined (0 0.5 0 0, 1 0 0 1, 2 0 1 1)" << endl;

    if(min_dust_dens == 0 && max_dust_dens == 0)
    {
        min_dust_dens = 0.1;
        max_dust_dens = 1;
    }

    if(min_dust_dens == 0)
    {
        min_dust_dens = 0.001 * max_dust_dens;
    }

    if(min_dust_dens / max_dust_dens > 0.9)
        min_dust_dens = 0.9 * max_dust_dens;

    point_fields[2] << "set cbrange[" << log10(min_dust_dens) << ":" << log10(max_dust_dens) << "]" << endl;
    point_fields[2] << "set format cb \'%.02g\'" << endl;

    point_fields[2] << "splot  '-' w p ls 1,'-' with vectors as 2,'-' with vectors as 1" << endl;

    // 3 gas_temp
    point_fields[3] << point_header.str();
    point_fields[3] << "set palette defined (0 0.05 0 0, 0.4 1 0 0, 0.7 1 1 0, 1 1 1 0.5)" << endl;

    point_fields[3] << "set title \'3D gas temperature distribution (min: " << min_gas_temp
                    << "[K]; max: " << max_gas_temp << "[K])\' font \'Arial,12\'" << endl;
    point_fields[3] << "set cblabel \'temperature [K]\'" << endl;

    if(min_gas_temp == 0 && max_gas_temp == 0)
    {
        min_gas_temp = 0.1;
        max_gas_temp = 1;
    }

    if(min_gas_temp / max_gas_temp > 0.90)
        min_gas_temp = 0.9 * max_gas_temp;

    point_fields[3] << "set cbrange[" << float(min_gas_temp) << ":" << float(max_gas_temp) << "]" << endl;
    point_fields[3] << "set format cb \'%.03g\'" << endl;

    point_fields[3] << "splot  '-' w p ls 1,'-' with vectors as 2,'-' with vectors as 1" << endl;

    // 4 dust temp
    point_fields[4] << point_header.str();
    point_fields[4] << "set palette defined (0 0.05 0 0, 0.4 1 0 0, 0.7 1 1 0, 1 1 1 0.5)" << endl;

    point_fields[4] << "set title \'3D dust temperature distribution (min: " << min_dust_temp
                    << "[K]; max: " << max_dust_temp << "[K])\' font \'Arial,12\'" << endl;
    point_fields[4] << "set cblabel \'temperature [K]\'" << endl;

    if(min_dust_temp == 0 && max_dust_temp == 0)
    {
        min_dust_temp = 0.1;
        max_dust_temp = 1;
    }

    if(min_dust_temp / max_dust_temp > 0.9)
        min_dust_temp = 0.9 * max_dust_temp;

    point_fields[4] << "set cbrange[" << float(min_dust_temp) << ":" << float(max_dust_temp) << "]" << endl;
    point_fields[4] << "set format cb \'%.03g\'" << endl;

    point_fields[4] << "splot  '-' w p ls 1,'-' with vectors as 2,'-' with vectors as 1" << endl;

    // 5 rat
    point_fields[5] << point_header.str();
    point_fields[5] << "set palette defined (0 0.05 0 0, 0.4 1 0 0, 0.7 1 1 0, 1 1 1 0.5)" << endl;

    point_fields[5] << "set title \'3D aligned radii distribution (min ID: " << aalg_min
                    << "; max ID: " << aalg_max << ")\' font \'Arial,12\'" << endl;

    point_fields[5] << "set cblabel \'aligned radius ID\'" << endl;

    if(aalg_min == aalg_max)
        aalg_max = 1.01 * aalg_min;

    point_fields[5] << "set cbrange[" << aalg_min << ":" << aalg_max << "]" << endl;

    point_fields[5] << "set format cb \'%.03g\'" << endl;

    point_fields[5] << "splot  '-' w p ls 1,'-' with vectors as 2,'-' with vectors as 1" << endl;

    vec_header.str("");
    vec_header << "reset" << endl;
    vec_header << "#set terminal postscript" << endl;
    vec_header << "#set output \'\'" << endl;
    vec_header << "set ticslevel 0" << endl;
    vec_header << "set size ratio -1" << endl;
    vec_header << "set view 45,45" << endl;

    vec_header << "set xlabel \'x[m]\'" << endl;
    vec_header << "set ylabel \'y[m]\'" << endl;
    vec_header << "set zlabel \'z[m]\'" << endl;

    vec_header << "set xrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;
    vec_header << "set yrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;
    vec_header << "set zrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;

    vec_header << "set style arrow 1 nohead ls 1 lw 1 lc rgb 0x0000cc" << endl;
    vec_header << "set style arrow 2 nohead ls 1 lw 1 lc rgb 0x5500dd" << endl;
    vec_header << "set style arrow 3 ls 1 lw 1 lc palette" << endl;

    vec_header << "set grid" << endl;
    vec_header << "set nokey" << endl;

    // 0 mag
    vec_fields[0] << vec_header.str();
    vec_fields[0] << "set palette defined (0 1 0 0, 0.5 0.0 0.9 0,  0.75 0.0 0.9 1, 0.9 0 0.1 0.9)" << endl;

    if(min_mag == 0 && max_mag == 0)
    {
        min_mag = 1e-45;
        max_mag = 2e-45;
    }

    vec_fields[0] << "set title \'3D mag. field distribution (min:" << log10(min_mag)
                  << " log10([T]); max:" << log10(max_mag) << " log10([T])  \' font \'Arial,12\'" << endl;

    if(min_mag / max_mag > 0.9)
        min_mag = 0.9 * max_mag;

    vec_fields[0] << "set cbrange[" << log10(min_mag) << ":" << log10(max_mag) << "]" << endl;
    vec_fields[0] << "set format cb \'%.02g\'" << endl;
    vec_fields[0] << "splot  \'-\' with vectors as 3, \'-\' with vectors as 2, \'-\' "
                     "with vectors as 1"
                  << endl;

    // 1 vel
    vec_fields[1] << vec_header.str();
    vec_fields[1] << "set palette defined (0 1 0 0, 0.5 0.0 0.9 0,  0.75 0.0 0.9 1, 0.9 0 0.1 0.9)" << endl;

    if(min_vel == 0 && max_vel == 0)
    {
        min_vel = 1e-45;
        max_vel = 1e-45;
    }

    vec_fields[1] << "set title \'3D vel. field directions (min:" << log10(min_vel)
                  << " log10(m/s); max:" << log10(max_vel) << " log10(m/s)\' font \'Arial,12\'" << endl;

    if(min_vel / max_vel > 0.9)
        min_vel = 0.9 * max_vel;

    vec_fields[1] << "set cbrange[" << float(log10(min_vel)) << ":" << float(log10(max_vel)) << "]" << endl;
    vec_fields[1] << "set format cb \'%.03g\'" << endl;
    vec_fields[1] << "splot  \'-\' with vectors as 3, \'-\' with vectors as 2, \'-\' "
                     "with vectors as 1"
                  << endl;

    for(uint pos = 0; pos < 8; pos++)
    {
        // point_fields[pos] << "\ne" << endl;
    }

    for(uint pos = 0; pos < 2; pos++)
    {
        vec_fields[pos] << "\ne" << endl;
    }

    stringstream buffer;
    buffer.str("");

    for(uint i_r = 0; i_r <= N_r; i_r++)
    {
        for(uint i_z = 0; i_z <= N_z; i_z++)
        {
            double r = listR[i_r];
            double z;
            if(i_r == N_r)
                z = listZ[N_r - 1][i_z];
            else
                z = listZ[i_r][i_z];

            double Nstep = double(int(15 + 25.0 * double(i_r) / double(N_r) + 0.5));
            double tmp_dph = PIx2 / Nstep;

            for(uint i_ph = 0; i_ph < Nstep; i_ph++)
            {
                Vector3D p1(r, tmp_dph * i_ph, z);
                Vector3D p2(r, tmp_dph * (i_ph + 1), z);

                p1.cyl2cart();
                p2.cyl2cart();
                Vector3D dist = p2 - p1;

                buffer << float(p1.X()) << " " << float(p1.Y()) << " " << float(p1.Z()) << " "
                       << float(dist.X()) << " " << float(dist.Y()) << " " << float(dist.Z()) << endl;
            }
        }
    } /**/

    for(uint pos = 0; pos < 8; pos++)
    {
        point_fields[pos] << buffer.str() << "\ne" << endl;
    }

    for(uint pos = 0; pos < 2; pos++)
    {
        vec_fields[pos] << buffer.str() << "\ne" << endl;
    }

    buffer.str("");

    for(uint i_r = 0; i_r < N_r; i_r++)
    {
        double dPh = PIx2 / double(N_ph[i_r]);
        for(uint i_ph = 0; i_ph <= N_ph[i_r]; i_ph++)
        {
            double r = listR[i_r];
            double ph = i_ph * dPh;

            Vector3D p1(r, ph, -Zmax);
            Vector3D p2(r, ph, Zmax);

            p1.cyl2cart();
            p2.cyl2cart();
            Vector3D dist = p2 - p1;

            buffer << float(p1.X()) << " " << float(p1.Y()) << " " << float(p1.Z()) << " " << float(dist.X())
                   << " " << float(dist.Y()) << " " << float(dist.Z()) << endl;
        }
    }

    for(uint pos = 0; pos < 8; pos++)
    {
        point_fields[pos] << buffer.str() << "\ne" << endl;
    }

    for(uint pos = 0; pos < 2; pos++)
    {
        vec_fields[pos] << "\ne" << endl;
    }

    for(uint i_z = 0; i_z <= N_z; i_z++)
    {
        double dPh = PIx2 / double(N_ph[0]);
        for(uint i_ph = 0; i_ph <= N_ph[0]; i_ph++)
        {
            double zmin = listZ[0][i_z];
            double zmax = listZ[N_r - 1][i_z];
            double ph = i_ph * dPh;

            Vector3D p1(Rmin, ph, zmin);
            Vector3D p2(Rmax, ph, zmax);

            p1.cyl2cart();
            p2.cyl2cart();
            Vector3D dist = p2 - p1;

            buffer << float(p1.X()) << " " << float(p1.Y()) << " " << float(p1.Z()) << " " << float(dist.X())
                   << " " << float(dist.Y()) << " " << float(dist.Z()) << endl;
        }
    }

    for(uint pos = 0; pos < 8; pos++)
    {
        point_fields[pos] << buffer.str();
    }

    for(uint pos = 0; pos < 8; pos++)
        point_fields[pos].close();

    for(uint pos = 0; pos < 2; pos++)
        vec_fields[pos].close();

    cout << "- Writing of plot files                : done" << endl;

    return true;
}

bool CGridCylindrical::saveBinaryGridFile(string filename, ushort id, ushort data_size)
{
    if(max_cells == 0)
    {
        cout << ERROR_LINE << "Cannot save cylindrical grid file to:" << endl;
        cout << filename;
        cout << "Not enough cells available! " << endl;
        return false;
    }

    ofstream bin_writer(filename.c_str(), ios::out | ios::binary);

    if(bin_writer.fail())
    {
        cout << ERROR_LINE << "Cannot open cylindrical grid file:" << endl;
        cout << filename;
        return false;
    }

    bin_writer.write((char *)&id, 2);
    bin_writer.write((char *)&data_size, 2);

    if(dataID == GRID_ID_CYL)
    {
        for(uint i = 0; i < data_offset; i++)
        {
            ushort tmp_ids = data_ids[i];
            bin_writer.write((char *)&tmp_ids, 2);
        }
    }
    else
    {
        cout << ERROR_LINE << "Cannot save cylindrical grid file to:" << endl;
        cout << filename;
        cout << "A cylindrical grid requires an ID of " << GRID_ID_CYL << "!" << endl;
        return false;
    }

    bin_writer.write((char *)&Rmin, 8);
    bin_writer.write((char *)&Rmax, 8);
    bin_writer.write((char *)&Zmax, 8);
    bin_writer.write((char *)&N_r, 2);
    bin_writer.write((char *)&N_ph[0], 2);
    if(log_factorZ == -1)
    {
        // Write N_z - 2 to take empty cells out
        uint tmp_N_z = N_z - 2;
        bin_writer.write((char *)&tmp_N_z, 2);
    }
    else
        bin_writer.write((char *)&N_z, 2);
    bin_writer.write((char *)&log_factorR, 8);
    bin_writer.write((char *)&log_factorPh, 8);
    bin_writer.write((char *)&log_factorZ, 8);
    if(log_factorR == 0)
        for(uint i_r = 1; i_r < N_r; i_r++)
            bin_writer.write((char *)&listR[i_r], 8);
    if(log_factorPh == 0)
        for(uint i_ph = 1; i_ph < N_ph[0]; i_ph++)
            bin_writer.write((char *)&listPh[0][i_ph], 8);
    else if(log_factorPh == -1)
        for(uint i_r = 0; i_r < N_r; i_r++)
            bin_writer.write((char *)&N_ph[i_r], 2);
    if(log_factorZ == 0)
        for(uint i_th = 1; i_th < N_z; i_th++)
            bin_writer.write((char *)&listZ[0][i_th], 8);
    else if(log_factorZ == -1)
        for(uint i_r = 0; i_r < N_r; i_r++)
        {
            double dz = listZ[i_r][2] - listZ[i_r][1];
            bin_writer.write((char *)&dz, 8);
        }

    for(uint i_r = 0; i_r < N_r; i_r++)
    {
        // cout << "-> Writing cylindrical grid file: " << float(100.0 * double(i_r) / double(N_r)) << "      \r"
        //      << flush;

        for(uint i_ph = 0; i_ph < N_ph[i_r]; i_ph++)
        {
            for(uint i_z = 0; i_z < N_z; i_z++)
            {
                if(log_factorZ != -1 || (i_z < N_z - 1 && i_z > 0))
                {
                    for(uint i = 0; i < data_offset; i++)
                    {
                        double tmp_data = grid_cells[i_r][i_ph][i_z]->getData(i);
                        bin_writer.write((char *)&tmp_data, 8);
                    }
                }
            }
        }
    }

    for(uint i_z = 0; i_z < N_z; i_z++)
    {
        if(log_factorZ != -1 || (i_z < N_z - 1 && i_z > 0))
        {
            for(uint i = 0; i < data_offset; i++)
            {
                double tmp_data = center_cells[i_z]->getData(i);
                bin_writer.write((char *)&tmp_data, 8);
            }
        }
    }

    bin_writer.close();

    cout << CLR_LINE;
    cout << "- Writing cylindrical grid file : done" << endl;

    return true;
}

void CGridCylindrical::printParameters()
{
    if(max_cells == 0)
        cout << ERROR_LINE << "No cylindrical grid parameters available! " << endl;
    else
    {
        cout << CLR_LINE;
        cout << "Cylindrical grid parameters (ID: " << getDataID() << ", data len.: " << getDataOffset()
             << ", Nr: " << N_r;
        if(log_factorPh == -1)
        {
            uint min_phi = MAX_UINT;
            uint max_phi = 0;
            for(uint i_r = 0; i_r < N_r; i_r++)
            {
                if(N_ph[i_r] < min_phi)
                    min_phi = N_ph[i_r];
                if(N_ph[i_r] > max_phi)
                    max_phi = N_ph[i_r];
            }
            cout << ", Nph (min,max): [" << min_phi << ", " << max_phi << "]";
        }
        else
            cout << ", Nph: " << N_ph[0];

        if(log_factorZ == -1)
            cout << ", Nz: " << N_z - 2 << ")" << endl;
        else
            cout << ", Nz: " << N_z << ")" << endl;
        cout << SEP_LINE;

        cout << "- Number of cylindrical cells   : " << max_cells << endl;

        printPhysicalParameters();
        cout << SEP_LINE;
    }
}

bool CGridCylindrical::positionPhotonInGrid(photon_package * pp)
{
    uint dirID = pp->getDirectionID();
    uint rID = MAX_UINT, zID = MAX_UINT, phID = MAX_UINT;

    if(dirID < 6 && dirID > 1 && pp->getPositionCell() != 0)
    {
        cell_cyl * tmp_cell = (cell_cyl *)pp->getPositionCell();

        rID = tmp_cell->getRID();
        zID = tmp_cell->getZID();
        phID = tmp_cell->getPhID();

        // Update index of next cell
        switch(dirID)
        {
            case 2:
                if(zID == 0)
                    return false;
                zID--;
                break;

            case 3:
                zID++;
                if(zID >= N_z)
                    return false;
                break;

            case 4:
                if(phID == 0)
                    phID += N_ph[rID];
                phID--;
                break;

            case 5:
                phID++;
                if(phID >= N_ph[rID])
                    phID -= N_ph[rID];
                break;

            default:
                return false;
                break;
        }

        // Set next cell
        if(rID == MAX_UINT)
            pp->setPositionCell(center_cells[zID]);
        else
            pp->setPositionCell(grid_cells[rID][phID][zID]);
        return true;
    }

    uint i_r = 0, i_ph = 0, i_z = 0;
    Vector3D pos = pp->getPosition();
    double cy_r = sqrt( pos.X()*pos.X() + pos.Y()*pos.Y() );

    if(cy_r < Rmin)
    {
        i_z = CMathFunctions::biListIndexSearch(pos.Z(), listZ[0], N_z + 1);
        if(i_z == MAX_UINT)
            return false;

        pp->setPositionCell(center_cells[i_z]);
        return true;
    }

    if(pp->getPositionCell() != 0 && (dirID == 0 || dirID == 1))
    {
        cell_cyl * tmp_cell = (cell_cyl *)pp->getPositionCell();
        rID = tmp_cell->getRID();

        if(dirID == 0)
            i_r = rID - 1;
        else
            i_r = rID + 1;
    }
    else
        i_r = CMathFunctions::biListIndexSearch(cy_r, listR, N_r + 1);

    if(i_r >= N_r || i_r == MAX_UINT)
        return false;

    if(N_ph[i_r] > 1)
    {
        double phi_tmp = pos.getPhiCoord();
        i_ph = CMathFunctions::biListIndexSearch(phi_tmp, listPh[i_r], N_ph[i_r] + 1);
        if(i_ph == MAX_UINT)
            return false;
    }

    i_z = CMathFunctions::biListIndexSearch(pos.Z(), listZ[i_r], N_z + 1);
    if(i_z == MAX_UINT)
        return false;

    pp->setPositionCell(grid_cells[i_r][i_ph][i_z]);
    return true;
}

bool CGridCylindrical::goToNextCellBorder(photon_package * pp)
{
    cell_cyl * tmp_cell = (cell_cyl *)pp->getPositionCell();
    Vector3D p = pp->getPosition();
    Vector3D d = pp->getDirection();

    bool hit = false;
    double path_length = 1e300;

    uint rID = tmp_cell->getRID();
    uint zID = tmp_cell->getZID();

    uint dirID = MAX_UINT;

    if(rID == MAX_UINT)
    {
        double A = d.X() * d.X() + d.Y() * d.Y();
        if(A > 0)
        {
            // --- Radial cell borders ---
            double r2 = Rmin * (1 + MIN_LEN_STEP*EPS_DOUBLE);

            double B = p.X() * d.X() + p.Y() * d.Y();
            double C = p.X() * p.X() + p.Y() * p.Y() - r2 * r2;
            // dscr is always positive, we are inside the inner cell
            double dscr = B * B - A * C;

            dscr = sqrt(dscr);
            double length = (-B + dscr) / A;

            if(length != 0 && length < path_length)
            {
                path_length = length;
                hit = true;
                dirID = 1;
            }
        }

        // --- vertical cell borders ---
        double z1 = listZ[0][zID] - (abs(listZ[0][zID]) + 1) * MIN_LEN_STEP*EPS_DOUBLE;
        double z2 = listZ[0][zID + 1] + (abs(listZ[0][zID + 1]) + 1) * MIN_LEN_STEP*EPS_DOUBLE;

        Vector3D v_n1(0, 0, -1);
        Vector3D v_a1(0, 0, z1);

        double den1 = v_n1 * d;
        // if den1 < 0, photon moves in +z, but p.Z() >= z1
        if(den1 > 0)
        {
            double num = v_n1 * (p - v_a1);
            double length = -num / den1;

            if(length != 0 && length < path_length)
            {
                path_length = length;
                hit = true;
                dirID = 2;
            }
        }

        Vector3D v_n2(0, 0, 1);
        Vector3D v_a2(0, 0, z2);

        double den2 = v_n2 * d;
        // if den2 < 0, photon moves in -z, but p.Z() <= z2
        if(den2 > 0)
        {
            double num = v_n2 * (p - v_a2);
            double length = -num / den2;

            if(length != 0 && length < path_length)
            {
                path_length = length;
                hit = true;
                dirID = 3;
            }
        }
    }
    else
    {
        // --- Radial cell borders ---
        double A = d.X() * d.X() + d.Y() * d.Y();
        if(A > 0)
        {
            // inner/outer cell border is reduced/enlarged to ensure
            // a large enough step length
            double r1 = listR[rID] * (1 - MIN_LEN_STEP*EPS_DOUBLE);
            double r2 = listR[rID + 1] * (1 + MIN_LEN_STEP*EPS_DOUBLE);

            double B = p.X() * d.X() + p.Y() * d.Y();
            double B_sq = pow(B, 2);

            double C1 = p.X() * p.X() + p.Y() * p.Y() - r1 * r1;
            double C2 = p.X() * p.X() + p.Y() * p.Y() - r2 * r2;

            double dscr1 = B_sq - A * C1;
            // dscr2 is always >= 0
            double dscr2 = B_sq - A * C2;

            if(dscr1 > 0)
            {
                dscr1 = sqrt(dscr1);
                // "+"-solution is not needed for inner cells; only the "-"-solution can be correct
                double length = (-B - dscr1) / A;

                if(length > 0 && length < path_length)
                {
                    path_length = length;
                    hit = true;
                    dirID = 0;
                }
            }

            dscr2 = sqrt(dscr2);
            // "-"-solution is not needed for outer cells; only the "+"-solution can be correct
            double length = (-B + dscr2) / A;

            if(length != 0 && length < path_length)
            {
                path_length = length;
                hit = true;
                dirID = 1;
            }
        }

        // --- vertical cell borders ---
        // inner/outer cell border is reduced/enlarged to ensure
        // a large enough step length
        double z1 = listZ[rID][zID] - (abs(listZ[rID][zID]) + 1) * MIN_LEN_STEP*EPS_DOUBLE;
        double z2 = listZ[rID][zID + 1] + abs(listZ[rID][zID + 1] + 1) * MIN_LEN_STEP*EPS_DOUBLE;

        // Vector3D v_n1(0, 0, -1);
        // Vector3D v_a1(0, 0, z1);

        //double den1 = v_n1 * d;
        double den1 = -d.Z();
        // if den1 < 0, photon moves in +z, but p.Z() >= z1
        if(den1 > 0)
        {
            //double num = v_n1 * (p - v_a1);
            double num = z1 - p.Z();
            double length = -num / den1;

            if(length != 0 && length < path_length)
            {
                path_length = length;
                hit = true;
                dirID = 2;
            }
        }

        // Vector3D v_n2(0, 0, 1);
        // Vector3D v_a2(0, 0, z2);

        //double den2 = v_n2 * d;
        double den2 = d.Z();
        // if den2 < 0, photon moves in -z, but p.Z() <= z2
        if(den2 > 0)
        {
            //double num = v_n2 * (p - v_a2);
            double num = p.Z() - z2;
            double length = -num / den2;

            if(length != 0 && length < path_length)
            {
                path_length = length;
                hit = true;
                dirID = 3;
            }
        }

        // --- Phi cell borders ---
        if(N_ph[rID] > 1)
        {
            uint phID = tmp_cell->getPhID();

            double r = sqrt(p.sq_length());

            double ph1 = listPh[rID][phID] * (1 - MIN_LEN_STEP*EPS_DOUBLE) - MIN_LEN_STEP*EPS_DOUBLE;
            double ph2 = listPh[rID][phID + 1] * (1 + MIN_LEN_STEP*EPS_DOUBLE) + MIN_LEN_STEP*EPS_DOUBLE;

            double sin_ph1 = sin(ph1);
            double sin_ph2 = sin(ph2);

            double cos_ph1 = cos(ph1);
            double cos_ph2 = cos(ph2);

            Vector3D v_n1 = -Vector3D(-sin_ph1, cos_ph1, 0);
            Vector3D v_a1 = r * Vector3D(cos_ph1, sin_ph1, 0);

            double den1 = v_n1 * d;
            if(den1 != 0)
            {
                double num = v_n1 * (p - v_a1);
                double length = -num / den1;

                if(length > 0 && length < path_length)
                {
                    path_length = length;
                    hit = true;
                    dirID = 4;
                }
            }

            Vector3D v_n2 = Vector3D(-sin_ph2, cos_ph2, 0);
            Vector3D v_a2 = r * Vector3D(cos_ph2, sin_ph2, 0);

            double den2 = v_n2 * d;
            if(den2 != 0)
            {
                double num = v_n2 * (p - v_a2);
                double length = -num / den2;

                if(length > 0 && length < path_length)
                {
                    path_length = length;
                    hit = true;
                    dirID = 5;
                }
            }
        }
    }

    if(!hit)
    {
        cout << ERROR_LINE << "Wrong cell border!                                   " << endl;
        return false;
    }

    pp->setPosition(p + d * path_length);

    if(p == pp->getPosition())
    {
        cout << ERROR_LINE << "Could not transfer photon to the next cell border!   " << endl;
        return false;
    }

    pp->setTmpPathLength(path_length);
    pp->setDirectionID(dirID);
    return true;
}

bool CGridCylindrical::updateShortestDistance(photon_package * pp)
{
    /*Vector3D tmp_pos;
    //    double min_dist, tmp_dist[6];

    //    double loc_x_min, loc_x_max, loc_y_min, loc_y_max, loc_z_min, loc_z_max;
    bool found = false;

    cell_oc * tmp_cell_pos = (cell_oc *)pp->getPositionCell();

    tmp_pos = pp->getPosition();

    loc_x_min = tmp_cell_pos->x_min;
    loc_y_min = tmp_cell_pos->y_min;
    loc_z_min = tmp_cell_pos->z_min;

    loc_x_max = tmp_cell_pos->x_max;
    loc_y_max = tmp_cell_pos->y_max;
    loc_z_max = tmp_cell_pos->z_max;

    min_dist = 1E200;

    tmp_dist[0] = abs(loc_x_max - tmp_pos.X());
    tmp_dist[1] = abs(tmp_pos.X() - loc_x_min);

    tmp_dist[2] = abs(loc_y_max - tmp_pos.Y());
    tmp_dist[3] = abs(tmp_pos.Y() - loc_y_min);

    tmp_dist[4] = abs(loc_z_max - tmp_pos.Z());
    tmp_dist[5] = abs(tmp_pos.Z() - loc_z_min);

    for(int i = 0; i < 6; i++)
    {
            if(min_dist > tmp_dist[i])
            {
                    min_dist = tmp_dist[i];
                    found = true;
            }
    }

    pp->setShortestDistance(min_dist);*/
    return false;
}

Vector3D CGridCylindrical::getCenter(const cell_basic & cell) const
{
    Vector3D center;
    const cell_cyl * tmp_cell = (const cell_cyl *)&cell;

    if(tmp_cell->getRID() == MAX_UINT)
    {
        double z = listZ[0][tmp_cell->getZID()];
        double dz = listZ[0][tmp_cell->getZID() + 1] - z;
        return Vector3D(0, 0, z + 0.5 * dz);
    }

    double r = listR[tmp_cell->getRID()];
    double dr = listR[tmp_cell->getRID() + 1] - r;
    double ph = listPh[tmp_cell->getRID()][tmp_cell->getPhID()];
    double dph = listPh[tmp_cell->getRID()][tmp_cell->getPhID() + 1] - ph;
    double z = listZ[tmp_cell->getRID()][tmp_cell->getZID()];
    double dz = listZ[tmp_cell->getRID()][tmp_cell->getZID() + 1] - z;

    double sin_ph = sin(ph + 0.5 * dph);
    double cos_ph = cos(ph + 0.5 * dph);

    center = (r + 0.5 * dr) * Vector3D(cos_ph, sin_ph, 0);
    center.setZ(z + 0.5 * dz);

    return center;
}

bool CGridCylindrical::next(photon_package * pp)
{
    if(!positionPhotonInGrid(pp))
        return false;

    if(!goToNextCellBorder(pp))
        return false;

    return true;
}

bool CGridCylindrical::findStartingPoint(photon_package * pp)
{
    double path_length = 1e300;
    bool hit = false;
    Vector3D p = pp->getPosition();
    Vector3D d = pp->getDirection();

    if(isInside(p))
        return positionPhotonInGrid(pp);

    // need distance to inner cell, we are outside the grid
    double r1 = Rmax * (1 - MIN_LEN_STEP*EPS_DOUBLE);
    double z1 = Zmax * (1 - MIN_LEN_STEP*EPS_DOUBLE);


    // --- Radial cell border ---
    double A = d.X() * d.X() + d.Y() * d.Y();
    if(A > 0)
    {
        double B = p.X() * d.X() + p.Y() * d.Y();
        double C = p.X() * p.X() + p.Y() * p.Y() - r1 * r1;
        double dscr = B * B - A * C;

        if(dscr > 0)
        {
            dscr = sqrt(dscr);
            // "+"-solution is not needed for inner cells; only the "-"-solution can be correct
            double length = (-B - dscr) / A;

            if(length >= 0 && length < path_length)
                if(abs(p.Z() + d.Z() * length) < z1)
                {
                    path_length = length;
                    hit = true;
                }
        }
    }

    // --- vertical cell borders ---
    Vector3D v_n1(0, 0, -1);
    double den1 = v_n1 * d;
    if(den1 > 0)
    {
        Vector3D v_a1(0, 0, z1);
        double num = v_n1 * (p - v_a1);
        double length = -num / den1;

        if(length != 0 && length < path_length)
            if(pow(p.X() + d.X() * length, 2) + pow(p.Y() + d.Y() * length, 2) < r1 * r1)
            {
                path_length = length;
                hit = true;
            }
    }

    Vector3D v_n2(0, 0, 1);
    double den2 = v_n2 * d;
    if(den2 > 0)
    {
        Vector3D v_a2(0, 0, -z1);
        double num = v_n2 * (p - v_a2);
        double length = -num / den2;

        if(length != 0 && length < path_length)
            if(pow(p.X() + d.X() * length, 2) + pow(p.Y() + d.Y() * length, 2) < r1 * r1)
            {
                path_length = length;
                hit = true;
            }
    }

    if(!hit)
        return false;

    pp->setPosition(p + d * path_length);
    pp->setDirectionID(MAX_UINT);
    return positionPhotonInGrid(pp);
}

void CGridCylindrical::getLengths(uint bins, double & step_xy, double & off_xy)
{
    step_xy = 2 * Rmax / double(bins);
    off_xy = step_xy / 2.0;
}

bool CGridCylindrical::createCellList()
{
    if(max_cells == 0)
    {
        cout << ERROR_LINE << "Cylindrical grid contains no cells!" << endl;
        cout << "       Cell list cannot be created!" << endl;
        return false;
    }

    cell_list = new cell_basic *[max_cells];
    ulong pos_counter = 0;

    // cout << CLR_LINE;
    // cout << "-> Creating cell list    : 0 [%]           \r" << flush;

    for(uint i_r = 0; i_r < N_r; i_r++)
    {
        // cout << "-> Creating cell list     : " << 100.0 * float(i_r) / float(N_r) << " [%]        \r"
        //      << flush;

        for(uint i_ph = 0; i_ph < N_ph[i_r]; i_ph++)
        {
            for(uint i_z = 0; i_z < N_z; i_z++)
            {
                cell_list[pos_counter] = (cell_basic *)grid_cells[i_r][i_ph][i_z];
                pos_counter++;
            }
        }
    }

    for(uint i_z = 0; i_z < N_z; i_z++)
    {
        cell_list[pos_counter] = (cell_basic *)center_cells[i_z];
        pos_counter++;
    }

    // cout << CLR_LINE;
    // cout << "- Creating of cell list         : done          \n" << flush;
    return true;
}

double CGridCylindrical::getVolume(const cell_basic & cell) const
{
    const cell_cyl * cell_pos = (const cell_cyl *)&cell;

    double volume = 0;

    if(cell_pos->getRID() == MAX_UINT)
    {
        double z1 = listZ[0][cell_pos->getZID()];
        double z2 = listZ[0][cell_pos->getZID() + 1];
        volume = PI * Rmin * Rmin * (z2 - z1);
    }
    else
    {
        double r1 = listR[cell_pos->getRID()];
        double r2 = listR[cell_pos->getRID() + 1];
        double ph1 = listPh[cell_pos->getRID()][cell_pos->getPhID()];
        double ph2 = listPh[cell_pos->getRID()][cell_pos->getPhID() + 1];
        double z1 = listZ[cell_pos->getRID()][cell_pos->getZID()];
        double z2 = listZ[cell_pos->getRID()][cell_pos->getZID() + 1];

        volume = 0.5 * (ph2 - ph1) * (r2 * r2 - r1 * r1) * (z2 - z1);
    }

    return volume;
}

Vector3D CGridCylindrical::rotateToCenter(const photon_package & pp, Vector3D dir, bool inv, bool phi_only) const
{
    const cell_cyl * cell_pos = (const cell_cyl *)pp.getPositionCell();
    double phi = pp.getPosition().getPhiCoord();

    double phi_center = cell_pos->getRID() == MAX_UINT
                            ? 0
                            : 0.5 * (listPh[cell_pos->getRID()][cell_pos->getPhID()] +
                                        listPh[cell_pos->getRID()][cell_pos->getPhID() + 1]);
    dir.rot(Vector3D(0, 0, 1), inv ? phi - phi_center : phi_center - phi);

    return dir;
}

bool CGridCylindrical::saveBinaryGridFile(string filename)
{
    return saveBinaryGridFile(filename, GRID_ID_CYL, data_offset);
}

bool CGridCylindrical::loadGridFromBinaryFile(parameters & param)
{
    return loadGridFromBinaryFile(param, 0);
}

void CGridCylindrical::clear()
{
    line_counter = 0;
    char_counter = 0;
    cout << "Final cleanup                          : done" << endl;
}

bool CGridCylindrical::getPolarRTGridParameter(double max_len,
                             double pixel_width,
                             uint max_subpixel_lvl,
                             dlist & _listR,
                             uint & N_polar_r,
                             uint *& N_polar_ph)
{
    return CGridBasic::getPolarRTGridParameterWorker(max_len,
                                                     pixel_width,
                                                     max_subpixel_lvl,
                                                     _listR,
                                                     N_polar_r,
                                                     N_polar_ph,
                                                     N_r,
                                                     listR);
}

bool CGridCylindrical::isInside(const Vector3D & pos) const
{
    double sq_r = pos.X() * pos.X() + pos.Y() * pos.Y();
    if(sq_r > Rmax * Rmax)
        return false;

    if(pos.Z() < -Zmax)
        return false;

    if(pos.Z() > Zmax)
        return false;

    return true;
}

void CGridCylindrical::setRndPositionInCell(photon_package * pp, CRandomGenerator * rand_gen)
{
    Vector3D pos;
    cell_cyl * tmp_cell = (cell_cyl *)pp->getPositionCell();
    double r1, r2, ph1, ph2, z1, z2;

    double rnd_r = rand_gen->getRND();
    double rnd_ph = rand_gen->getRND();
    double rnd_z = rand_gen->getRND();

    if(tmp_cell->getZID() == MAX_UINT)
        return;
    else if(tmp_cell->getRID() == MAX_UINT)
    {
        r1 = 0;
        r2 = listR[0];
        ph1 = listPh[0][0];
        ph2 = listPh[0][N_ph[0]];
        z1 = listZ[0][0];
        z2 = listZ[0][N_z];
    }
    else
    {
        r1 = listR[tmp_cell->getRID()];
        r2 = listR[tmp_cell->getRID() + 1];
        ph1 = listPh[tmp_cell->getRID()][tmp_cell->getPhID()];
        ph2 = listPh[tmp_cell->getRID()][tmp_cell->getPhID() + 1];
        z1 = listZ[tmp_cell->getRID()][tmp_cell->getZID()];
        z2 = listZ[tmp_cell->getRID()][tmp_cell->getZID() + 1];
    }

    double sin_ph = sin(ph1 + rnd_ph * (ph2 - ph1));
    double cos_ph = cos(ph1 + rnd_ph * (ph2 - ph1));

    pos = (pow(pow(r1, 2) + rnd_r * (pow(r2, 2) - pow(r1, 2)), 1.0 / 2.0) * Vector3D(cos_ph, sin_ph, 0)) +
            Vector3D(0, 0, z1 + rnd_z * (z2 - z1));

    pp->setPosition(pos);
}
