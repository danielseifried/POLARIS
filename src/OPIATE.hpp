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
/* History:   05.12.2016                                                    */
/****************************************************************************/

#pragma once

#include <CCfits/CCfits>

#include "CommandParser.hpp"
#include "GridBasic.hpp"
#include "Parameters.hpp"


using namespace std;
using namespace CCfits;

#ifndef OPIATE_H
#define OPIATE_H

class COpiateDataBase
{
  public:
    COpiateDataBase()
    {
        current_index = MAX_UINT;
        list_freq = 0;
        list_IDs = 0;
        list_weight = 0;
        max_species = 0;
        max_ids = 0;

        database_counter = 0;

        max_velocity = 1;
        velocity_channel = 0;
        
        path_emi="";
        path_abs="";
        
        has_emi_data=false;
        has_abs_data=false;
    }

    ~COpiateDataBase()
    {
        cout << CLR_LINE;
        cout << "Final OPIATE data cleanup ...   \r" << flush;

        if(list_freq != 0)
        {
            delete[] list_freq;
            list_freq = 0;
        }

        if(list_IDs != 0)
        {
            delete[] list_IDs;
            list_IDs = 0;
        }

        if(list_weight != 0)
        {
            delete[] list_weight;
            list_weight = 0;
        }

        if(velocity_channel != 0)
        {
            delete[] velocity_channel;
            velocity_channel = 0;
        }

        if(entries.size()>0)
        {
            for(uint i =0;i<entries.size();i++)
            {
                COpiateEntry * tmp=entries[i];
                delete tmp;
                tmp=0;
                entries[i]=0;
            }
            
            entries.clear();
        }

        cout << CLR_LINE;
    }

    double getGaussA(double temp_gas, double v_turb);
    
    void calcLineBroadening(CGridBasic * grid);
    
    double getProjCellVelocityInterp(const Vector3D & tmp_pos,
                                     const Vector3D & dir_map_xyz,
                                     const VelFieldInterp & vel_field_interp);
    
    uint getIndex(uint op_id) const;

    void getMatrices(CGridBasic * grid,
                                 const photon_package * pp,
                                 uint i_spec,
                                 uint i_trans,
                                 double velocity,
                                 const LineBroadening & line_broadening,
                                 const MagFieldInfo & mag_field_info,
                                 StokesVector * line_emissivity,
                                 Matrix2D * line_absorption_matrix) const;
    
    uint getMaxSpecies();

    double getFrequency(uint pos);

    bool readDataBase(string filename);

    bool readOpiateDataBase(parameters & param);

    void printParameters(parameters & param, CGridBasic * grid);

    bool findIndexByName(string name);
    
    string getCurrentName();

    bool readEmissivityData(string filename);

    bool readAbsorptionData(string filename);

    bool readFitsData(string filename, Matrix2D & mat);

    double getCurrentFrequency();

    double getMolecularWeight();

    bool initVelChannels(uint nr_of_channels, double max_vel);

    double getVelocityChannel(uint vch);

    uint getNrOfVelChannels();

    double getMaxVelocity();

  private:
    class COpiateEntry
    {
      public:
        COpiateEntry()
        {
            freq=0;
            weight=0;
            size=0;
            name="Empty";
            em=0;
            ex=0;
        }

        COpiateEntry(uint _size)
        {
            freq=0;
            weight=0;
            size=_size;
            name="Empty";
            em=new double[size];
            ex=new double[size];

            for(uint i=0;i<size;i++)
            {
                em[i]=0;
                ex[i]=0;
            }               
        }

        ~COpiateEntry()
        {
            if(em != 0)
            {
                delete[] em;
                em = 0;
            }

            if(ex != 0)
            {
                delete[] ex;
                ex = 0;
            }
        }

        void setData(uint pos, double _em, double _ex)
        {
            em[pos]=_em;
            ex[pos]=_ex;
        }

        double freq;
        double weight;
        uint size;
        string name;
        double * em;
        double * ex;
    };    

    Matrix2D mat_emissivity;
    Matrix2D mat_absorption;

    uint current_index;
    uint max_species;
    uint max_ids;
    
    string path_emi;
    string path_abs;

    double * list_freq;
    double * list_IDs;
    double * list_weight;
    strlist list_names;

    double * velocity_channel;
    uint nr_of_velocity_channels;
    double max_velocity;

    uint database_counter;
    
    bool has_emi_data;
    bool has_abs_data;
    
    vector < COpiateEntry * > entries;
};

#endif
