#include <iostream>
#include <fstream>
#include <map>
#include <utility>
#include <string>
#include <sstream>
#include <exception>
#include <tuple>
#include <array>
#include <cmath>
#include <iomanip>
#include <boost/filesystem.hpp>
#include <boost/range/iterator_range.hpp>


constexpr char pesheader[] = "#          beta[]       gamma[deg]      energy[MeV]"
                             "   Bx[hbar^2/MeV]   By[hbar^2/MeV]   Bz[hbar^2/MeV]"
                             "  Bbb[hbar^2/MeV]  Bbg[hbar^2/MeV]  Bgg[hbar^2/MeV]"
                             "        Vrot[MeV]        Vvib[MeV]       Vcoll[MeV]";

using PEStable = std::map< std::pair<double,double> , std::array<double,10> >;

template <class ISTREAM>
PEStable readPES( ISTREAM & in )
{
    PEStable pestable;

    std::string line;
    while( std::getline(in,line) )
    {
        if( line == pesheader ) continue;
        
        double beta;
        double gamma;
        double energy;
        double Bx;
        double By;
        double Bz;
        double Bbb;
        double Bbg;
        double Bgg;
        double Vrot;
        double Vvib;
        double Vcoll;

        std::stringstream(line) >> beta >> gamma >> energy
                                >>   Bx >>    By >>     Bz
                                >>  Bbb >>   Bbg >>    Bgg
                                >> Vrot >>  Vvib >> Vcoll;
        
        std::pair<double,double> beta_gamma = std::make_pair( beta , gamma );
        std::array<double,10>    values     = { energy , Bx , By , Bz , Bbb , Bbg , Bgg , Vrot , Vvib , Vcoll };
        
        pestable[ beta_gamma ] = values;
    }

    return pestable;
}

template <class OSTREAM>
void writePES( PEStable const & pestable , OSTREAM & out )
{
    out << pesheader << std::endl;

    for( auto const & [ beta_gamma , values ] : pestable )
    {
        auto const [beta,gamma]                                  = beta_gamma;
        auto const [energy,Bx,By,Bz,Bbb,Bbg,Bgg,Vrot,Vvib,Vcoll] = values;


        out.precision(3);
        out << std::setw(17) << std::showpos << std::fixed;
        out << beta;
        
        out.precision(3);
        out << std::setw(17) << std::showpos << std::fixed;
        out << gamma;

        out.precision(4);
        out << std::setw(17) << std::showpos << std::fixed;
        out << energy;

        out.precision(4);
        out << std::setw(17) << std::showpos << std::fixed;
        out << Bx;
        
        out.precision(4);
        out << std::setw(17) << std::showpos << std::fixed;
        out << By;

        out.precision(4);
        out << std::setw(17) << std::showpos << std::fixed;
        out << Bz;

        out.precision(4);
        out << std::setw(17) << std::showpos << std::fixed;
        out << Bbb;
        
        out.precision(4);
        out << std::setw(17) << std::showpos << std::fixed;
        out << Bbg;

        out.precision(4);
        out << std::setw(17) << std::showpos << std::fixed;
        out << Bgg;

        out.precision(4);
        out << std::setw(17) << std::showpos << std::fixed;
        out << Vrot;

        out.precision(4);
        out << std::setw(17) << std::showpos << std::fixed;
        out << Vvib;

        out.precision(4);
        out << std::setw(17) << std::showpos << std::fixed;
        out << Vcoll;


        out << std::endl;
    }

    return;
}

template <class ISTREAM>
std::pair< std::pair<double,double> , std::array<double,10> > read_dirhbout( ISTREAM & in )
{
    double beta;
    double gamma;
    double energy;
    double Bx;
    double By;
    double Bz;
    double Bbb;
    double Bbg;
    double Bgg;
    double Vrot;
    double Vvib;

    bool beta_flag   = false;
    bool gamma_flag  = false;
    bool energy_flag = false;
    bool Bx_flag     = false;
    bool By_flag     = false;
    bool Bz_flag     = false;
    bool Bbb_flag    = false;
    bool Bbg_flag    = false;
    bool Bgg_flag    = false;
    bool Vrot_flag   = false;
    bool Vvib_flag   = false;



    std::string line;
    while( std::getline(in,line) )
    {
        if( line.rfind( " beta ................"                , 0 ) == 0 )
        {
            double betan;
            double betap;
            line.erase(0, std::string(" beta ................").size() );
            std::stringstream(line) >> betan >> betap >> beta;
            beta_flag = true;
            continue;
        }
        if( line.rfind( " gamma................"                , 0 ) == 0 )
        {
            double gamman;
            double gammap;
            line.erase(0, std::string(" gamma................").size() );
            std::stringstream(line) >> gamman >> gammap >> gamma;
            gamma_flag = true;
            continue;
        }
        if( line.rfind( " Total Energy ........"                , 0 ) == 0 )
        {
            line.erase(0, std::string(" Total Energy ........").size() );
            std::stringstream(line) >> energy;
            energy_flag = true;
            continue;
        }
        if( line.rfind( "...................Bx................" , 0 ) == 0 )
        {
            line.erase(0, std::string("...................Bx................").size() );
            std::stringstream(line) >> Bx;
            Bx_flag = true;
            continue;
        }
        if( line.rfind( "...................By................" , 0 ) == 0 )
        {
            line.erase(0, std::string("...................By................").size() );
            std::stringstream(line) >> By;
            By_flag = true;
            continue;
        }
        if( line.rfind( "...................Bz................" , 0 ) == 0 )
        {
            line.erase(0, std::string("...................Bz................").size() );
            std::stringstream(line) >> Bz;
            Bz_flag = true;
            continue;
        }
        if( line.rfind( " Mass parameters (Bbb,Bbg,Bgg)......." , 0 ) == 0 )
        {
            line.erase(0, std::string(" Mass parameters (Bbb,Bbg,Bgg).......").size() );
            std::stringstream(line) >> Bbb >> Bbg >> Bgg;
            Bbb_flag = true;
            Bbg_flag = true;
            Bgg_flag = true;
            continue;
        }
        if( line.rfind( " Rotational correction..............." , 0 ) == 0 )
        {
            line.erase(0, std::string(" Rotational correction...............").size() );
            std::stringstream(line) >> Vrot;
            Vrot_flag = true;
            continue;
        }
        if( line.rfind( " Vibrational correction.............." , 0 ) == 0 )
        {
            line.erase(0, std::string(" Vibrational correction..............").size() );
            std::stringstream(line) >> Vvib;
            Vvib_flag = true;
            continue;
        }
    }



    if( !beta_flag || !gamma_flag || !energy_flag ||
        !Bx_flag   || !By_flag    || !Bz_flag     ||
        !Bbb_flag  || !Bbg_flag   || !Bgg_flag    ||
        !Vrot_flag || !Vvib_flag 
      )
        throw std::runtime_error("{beta,gamma,energy,Bx,By,Bz,Bbb,Bbg,Bgg,Vrot,Vvib,Vcoll} not found!");


    std::pair<double,double> beta_gamma = std::make_pair( beta , gamma );
    std::array<double,10>    values     = { energy , Bx , By , Bz , Bbb , Bbg , Bgg , Vrot , Vvib , energy-Vrot-Vvib };

    return std::make_pair( beta_gamma , values );
}

template <class ISTREAM>
bool converged_dirhbout( ISTREAM & in )
{

    std::string line;
    while( std::getline(in,line) )
    {
        if( line.find("interrupted") != std::string::npos )
            return false;

        if( line.find("converged") != std::string::npos )
            return true;
    }

    throw std::runtime_error("neither converged nor interrupted!");
}


int main(int argc, char *argv[])
{
    try
    {
        /** Read existing pes.out data **/
        std::fstream pesout;
        pesout.open( "../../pes.out" , std::fstream::in | std::fstream::out |  std::fstream::app );
        if( !pesout )
        {
            pesout.open( "../../pes.out" ,  std::fstream::in | std::fstream::out );
            if( pesout.is_open() == false )
            {
                std::cerr << "pes.out opening failed!" << std::endl;
                return -1;
            }
        }

        auto pestable = readPES(pesout);
        pesout.close();



        /** Read new data **/
        if( boost::filesystem::is_directory("../../output") )
        {
            for( auto const & entry : boost::make_iterator_range( boost::filesystem::directory_iterator("../../output"), {} ) )
            {
                std::string folder = entry.path().string();

                std::ifstream dirhbout( folder + "/dirhb.out" );

                if( dirhbout.is_open() == false )
                    std::cerr << "Folder: " << folder << ", doesn't contain dirhb.out file!" << std::endl;
                else
                {
                    try
                    {
                        if( converged_dirhbout( dirhbout ) == false )
                            throw std::runtime_error("convergence failed!");
                        
                        // rewind dirhbout
                        dirhbout.clear();
                        dirhbout.seekg(0);

                        auto const [ beta_gamma , values ] = read_dirhbout( dirhbout );
                    
                        dirhbout.close();


                        // rounding to three digits after decimal point to avoid rounding-error
                        // troubles when insertion and comparison inside a map
                        auto [ beta , gamma ] = beta_gamma;
                        beta  = std::round( beta  * 1000.0 ) / 1000.0;
                        gamma = std::round( gamma * 1000.0 ) / 1000.0;


                        pestable[ std::pair<double,double>(beta,gamma) ] = values;
                    }
                    catch( std::exception const & e )
                    {
                        std::cout << "Exception in " << folder + "/dirhb.out: " << e.what() << std::endl;
                    }

                }

            }
        }



        /** Dump to pes.out file **/
        pesout.open( "../../pes.out" , std::fstream::out | std::fstream::trunc );
        writePES( pestable , pesout );
        pesout.close();

    }
    catch( std::exception const & e )
    {
        std::cout << "Exception: " << e.what() << std::endl;
    }

    return 0;
}
