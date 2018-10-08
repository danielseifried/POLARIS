#pragma once
#include "typedefs.h"
#include "chelper.h"
#include "Source.h"
#include "Detector.h"
#include "Dust.h"
#include "GasSpecies.h"
#include "CommandParser.h"

#ifndef CPIPELINE
#define CPIPELINE

class CPipeline
{
public:

    CPipeline(void)
    {
        begin = 0;
        end = 0;
        len = 0;

        m = 0;
        s = 0;
        h = 0;

        ru[0] = '|';
        ru[1] = '/';
        ru[2] = '-';
        ru[3] = '\\';
    }

    ~CPipeline(void)
    {
    }

    bool calcMonteCarloRadiationField(parameter & param);
    bool calcPolarizationMapsViaRayTracing(parameter & param);
    bool calcPolarizationMapsViaSynchrotron(parameter & param);
    bool calcChMapsViaRayTracing(parameter & param);
    bool calcPolarizationMapsViaMC(parameter & param);
    bool proberobeLineOfSight(parameter & param);

    bool assignDustMixture(parameter & param, CDustMixture * dust, CGridBasic * grid);
    bool assignGasSpecies(parameter & param, CGasMixture * gas, CGridBasic * grid);

    //bool calcRadPressure(parameter & param);

    bool solveNextVelocityField(uint itID)
    {
        return true;
    };
    //bool preparePressureData(CGridBasic * grid, CDustMixture * dust, parameter & param, bool plot, uint itID);

    bool createWavelengthList(parameter & param, CDustMixture * dust, CGasMixture * gas = 0);

    void printParameter(parameter & param, uint max_id);

    bool Init(int argc, char** argv);
    void Finish();
    void Error();
    void Run();

    bool createOutputPaths(string path)
    {
        const char * sep = SEP;
        string::size_type pos1 = 0, pos2 = 0;
        int len = -1;
        uint offset;

#ifdef WINDOWS
        offset = 0;
#else
        offset = 1;
#endif

        path_data = path + "data" + sep;
        path_plot = path + "plots" + sep;

        strlist folders;
        string folder;

        while(path.find(sep, offset) != string::npos)
        {
            pos1 = path.find(sep, offset);
            pos2 = path.find(sep, pos1 + 1);

            len = int(pos2 - pos1 + 1);

            if(len < 0) break;

            folder = path.substr(pos1 + 1, len - 2);
            path.erase(pos1 + 1, len - 1);
            folders.push_back(folder);
        }

        for(uint i = 0; i < folders.size(); i++)
        {
            folder = folders[i];
            path += folder + sep;

            if(!createPath(path))
            {
                cout << "ERROR: Failed to create output folder for data!" << endl;
                cout << path << std::endl;
                return false;
            }
        }

        if(!createPath(path))
        {
            cout << "ERROR: Failed to create output folder(s)!" << endl;
            cout << path << std::endl;
            return false;
        }

        if(!createPath(path_data))
        {
            cout << "ERROR: Failed to create output folder for data!" << endl;
            cout << path_data << std::endl;
            return false;
        }

        if(!createPath(path_plot))
        {
            cout << "ERROR: Failed to create output folder for plots!" << endl;
            cout << path_plot << std::endl;
            return false;
        }

        return true;
    }

    bool assignGridType(CGridBasic * &grid, parameter & param);

    bool createPath(string path)
    {
        const char * tmp_path = path.c_str();
#ifdef WINDOWS
        if(_access(tmp_path, 0) != 0)
            if(_mkdir(tmp_path)) return false;
#else
        if(access(tmp_path, 0) != 0)
            if(mkdir(tmp_path, S_IRWXU | S_IRWXG | S_IRWXO))
                return false;
#endif

        return true;
    }

    void printConversionParameter(parameter & param)
    {
        cout << "Conversion factors" << endl;
        cout << "- Conv. gas dens. in SI        : " << param.getSIConvDH() << endl;
        cout << "- Conv. length in SI           : " << param.getSIConvLength() << endl;
        cout << "- Conv. vel. field in SI.      : " << param.getSIConvVField() << endl;
        cout << "- Conv. mag. field in SI.      : " << param.getSIConvBField() << endl;
        cout << "- Mass fraction (Mdust/Mgas)   : " << param.getMassFraction() << endl;
        cout << "- Relative molecular mass (mu) : " << param.getMu() << endl;
    }

    void printPathParameter(parameter & param)
    {
        cout << "- Path grid file : " << param.getPathGrid() << endl;
        cout << "- Path output    : " << param.getPathOutput() << endl;
        cout << "- Nr. of threads : " << param.getNrOfThreads() << endl;
    }

    void printPlotParameter(parameter & param, bool input_output=false)
    {
        if(param.getNrOfGnuPoints() + param.getNrOfGnuVectors() + param.getInpMidDataPoints() +
                param.getOutMidDataPoints() + param.getInpAMIRAPoints() + param.getOutAMIRAPoints() > 0)
            cout << "Plot parameter" << endl;

        if(param.getNrOfGnuPoints() + param.getNrOfGnuVectors() != 0)
        {
            cout << "- Gnuplot                      : " << param.getNrOfGnuPoints() << " points, ";
            cout << param.getNrOfGnuVectors() << " vectors, ";
            cout << param.getmaxGridLines() << " max. level" << endl;
        }

        if(param.getInpMidDataPoints() != 0 && param.getOutMidDataPoints() != 0)
        {
            cout << "- Midplane            (in,out) : " << param.getInpMidDataPoints()
                    << ", " << param.getOutMidDataPoints() << " pixel" << endl;
        }
        else
        {
            if(param.getInpMidDataPoints() != 0)
                cout << "- Midplane                (in) : " << param.getInpMidDataPoints() << " pixel" << endl;

            if(param.getOutMidDataPoints() != 0)
                cout << "- Midplane               (out) : " << param.getOutMidDataPoints() << " pixel" << endl;
        }

        if(param.getInpAMIRAPoints() != 0 && param.getOutAMIRAPoints() != 0)
            cout << "- Amira points (in,out)        : " << param.getInpAMIRAPoints()
                << ", " << param.getOutAMIRAPoints() << " pixel" << endl;
        else
        {
            if(param.getInpAMIRAPoints() != 0)
                cout << "- Amira points            (in) : " << param.getInpAMIRAPoints() << " pixel" << endl;
            if(param.getOutAMIRAPoints() != 0)
                cout << "- Amira points           (out) : " << param.getOutAMIRAPoints() << " pixel" << endl;
        }
    }

    void printAlignmentParameter(parameter & param)
    {
        cout << "Dust grain alignment" << endl;

        if(param.getAligRANDOM())
            cout << "- Random alignment" << endl;
        else if(param.getAligPA())
            cout << "- Perfect alignment" << endl;
        else
        {
            if(param.getAligIDG())
                cout << "- IDG (paramagnetic alignment)" << endl;
            if(param.getAligRAT())
                cout << "- RAT (radiative torques); f_highJ: " << param.getFHighJ() << endl;
            if(param.getAligGOLD())
                cout << "- GOLD (mechanical alignment)" << endl;
            if(param.getAligINTERNAL())
                cout << "- Internal disalignment; f_c: " << param.getFcorr() << endl;
        }
    }

    void printSourceParameter(parameter & param, bool only_background=false, bool no_dust=false)
    {
        if(param.getNrOfSources() > 0)
        {
            if(param.getNrOfSources() != param.getNrOfPointSources() ||
                (param.getNrOfSources() == param.getNrOfPointSources() && 
                (!only_background || param.getScatteringToRay())))
                cout << "Defined radiation sources" << endl;

            if(param.getNrOfPointSources() > 0 &&
                    (!only_background || param.getScatteringToRay()))
                cout << "- Star(s)        : " << param.getNrOfPointSources() << endl;

            if(param.getNrOfDiffuseSources() > 0 && !only_background)
                cout << "- Star field(s)  : " << param.getNrOfDiffuseSources() << endl;

            if(param.getNrOfBackgroundSources() > 0)
                cout << "- Background(s)  : " << param.getNrOfBackgroundSources() << endl;

            if(param.getDustSource() && !only_background && !no_dust)
                cout << "- Dust as sources: yes" << endl;

            if(param.getISRFSource() && !only_background)
                cout << "- ISRF as sources: yes" << endl;
        }
    }

    void printDetectorParameter(parameter & param, bool monte_carlo=false)
    {
        if(monte_carlo)
        {
            cout << "Monte-Carlo parameter" << endl;
            if(param.getNrOfDustMCDetectors() > 0)
                cout << "- Nr. of detectors      : " << param.getNrOfDustMCDetectors() << endl;
        }
        else
        {
            cout << "Raytrace parameter" << endl;
            if(param.getNrOfDustRayDetectors() > 0)
                cout << "- Nr. of detectors      : " << param.getNrOfDustRayDetectors() << endl;
            param.printRTGridDescription();
            cout << "- Start, Stop           : " << param.getStart() + 1
                << ", " << param.getStop() + 1 << endl;
        }
        cout << "- Observer distance [m] : " << param.getMinObserverDistance()
            << " - " << param.getMaxObserverDistance() << endl;
        cout << "- Pixel in x-dir/nsides : " << param.getMinDetectorPixelX() << " (min) - "
            << param.getMaxDetectorPixelX() << " (max)" << endl;
        if(param.getMinDetectorPixelY() != 0 || param.getMaxDetectorPixelY() != 0)
            cout << "- Pixel in y-direction  : " << param.getMinDetectorPixelY() << " (min) - "
                << param.getMaxDetectorPixelY() << " (max)" << endl;
        if(!monte_carlo)
            cout << "- Max subpixel level    : " << param.getMaxSubpixelLvl() << endl;

        if(param.getMaxSidelengthX() > 0)
        {
            cout << "- Sidelength X [m]      : " << param.getMinSidelengthX() << " (min) - "
                    << param.getMaxSidelengthX() << " (max)" << endl;
            if(param.getUseGridSidelengthX())
                cout << "  (some detectors are using the grid size)" << endl;
        }
        else
            if(param.getUseGridSidelengthX())
                cout << "- Sidelength X [m]      : only the grid size is used" << endl;
        if(param.getMaxSidelengthY() > 0)
        {
            cout << "- Sidelength Y [m]      : " << param.getMinSidelengthY() << " (min) - "
                    << param.getMaxSidelengthY() << " (max)" << endl;
            if(param.getUseGridSidelengthY())
                cout << "  (some detectors are using the grid size)" << endl;
        }
        else
            if(param.getUseGridSidelengthY())
                cout << "- Sidelength Y [m]      : only the grid size is used" << endl;

        if(!monte_carlo)
        {
            if(param.getMaxMapShiftX() > 0)
                cout << "- Map shift X [m]       : " << param.getMinMapShiftX() << " (min) - "
                        << param.getMaxMapShiftX() << " (max)" << endl;
            if(param.getMaxMapShiftY() > 0)
                cout << "- Map shift Y [m]       : " << param.getMinMapShiftY() << " (min) - "
                        << param.getMaxMapShiftY() << " (max)" << endl;
        }

        cout << "- Rotation axis         : n1 = ";
        param.getAxis1().printValues();
        cout << ", angle [°] = " << param.getMinDetectorAngle1() << " (min) - "
                << param.getMaxDetectorAngle1() << " (max)" << endl;
        cout << "- Rotation axis         : n2 = ";
        param.getAxis2().printValues();
        cout << ", angle [°] = " << param.getMinDetectorAngle2() << " (min) - "
                << param.getMaxDetectorAngle2() << " (max)" << endl;
    }

    void printSynchrotronParameter(parameter & param)
    {
        cout << SEP_LINE;
        cout << "Synchotron parameter" << endl;
        cout << SEP_LINE;
        cout << "Observed wavelengths" << endl;
        dlist sync_ray_detectors = param.getSyncRayDetectors();
        for(uint i = 0; i < sync_ray_detectors.size(); i += NR_OF_RAY_DET)
        {
            uint pos = i / NR_OF_RAY_DET;
            double lam_min = sync_ray_detectors[i + 0];
            double lam_max = sync_ray_detectors[i + 1];
            uint lam_skip = uint(sync_ray_detectors[i + 2]);
            cout << "    - Synchotron emission detector " << (pos + 1) << ": from wl = " << lam_min << " [m] to wl = "
                << lam_max << " [m] with " << lam_skip << " steps" << endl;
        }
    }

    void deleteSourceLists()
    {
        for(uint s = 0; s < sources_mc.size(); s++)
            delete sources_mc[s];
        sources_mc.clear();

        for(uint s = 0; s < sources_ray.size(); s++)
            delete sources_ray[s];
        sources_ray.clear();
    }

private:
    bool writeSources(parameter & param, CGridBasic * grid);
    void createSourceLists(parameter & param, CDustMixture * dust, CGridBasic * grid);
    CDetector * createDetectorList(parameter & param, CDustMixture * dust, CGridBasic * grid);
    string path_plot, path_data;
    slist sources_mc, sources_ray;
    long h, m, s;
    double begin, end;

    long len;
    parameter_list param_list;
    unsigned char ru[4];

    //uint gria;
};

#endif
