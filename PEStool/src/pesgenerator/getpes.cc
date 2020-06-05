#include <iostream>
#include <fstream>
#include <map>
#include <utility>
#include <string>
#include <sstream>
#include <exception>
#include <tuple>
#include <cmath>
#include <iomanip>
#include <boost/filesystem.hpp>
#include <boost/range/iterator_range.hpp>


constexpr char pesheader[] = "#             beta[]          gamma[deg]         energy[MeV]";

using PEStable = std::map< std::pair<double,double> , double >;

template <class ISTREAM>
PEStable readPES( ISTREAM & in )
{
    PEStable pestable;

    std::string line;
    while( std::getline(in,line) )
    {
        if( line == pesheader ) continue;
        double beta, gamma, energy;
        std::stringstream(line) >> beta >> gamma >> energy;
        pestable[ std::pair<double,double>(beta,gamma) ] = energy;
    }

    return pestable;
}

template <class OSTREAM>
void writePES( PEStable const & pestable , OSTREAM & out )
{
    out << pesheader << std::endl;

    for( auto const & [ beta_gamma , energy ] : pestable )
    {
        double beta  = beta_gamma.first;
        double gamma = beta_gamma.second;

        out.precision(3);
        out << std::setw(20) << std::showpos << std::fixed;
        out << beta;
        out << std::setw(20) << std::showpos << std::fixed;
        out << gamma;
        out.precision(4);
        out << std::setw(20) << std::showpos << std::fixed;
        out << energy;
        out << std::endl;
    }

    return;
}


template <class ISTREAM>
std::tuple<double,double,double> read_dirhbout( ISTREAM & in )
{
    double beta;
    double gamma;
    double energy;
    
    bool   betaflag = false;
    bool  gammaflag = false;
    bool energyflag = false;
    
    std::string line;
    while( std::getline(in,line) )
    {
        if( line.rfind( " beta ................" , 0 ) == 0 )
        {
            double betan;
            double betap;
            line.erase(0, std::string(" beta ................").size() );
            std::stringstream(line) >> betan >> betap >> beta;
            betaflag = true;
            continue;
        }
        if( line.rfind( " gamma................" , 0 ) == 0 )
        {
            double gamman;
            double gammap;
            line.erase(0, std::string(" gamma................").size() );
            std::stringstream(line) >> gamman >> gammap >> gamma;
            gammaflag = true;
            continue;
        }
        if( line.rfind( " Total Energy ........" , 0 ) == 0 )
        {
            line.erase(0, std::string(" Total Energy ........").size() );
            std::stringstream(line) >> energy;
            energyflag = true;
            continue;
        }
    }
    
    if( betaflag==false || gammaflag==false || energyflag==false )
        throw std::runtime_error("Error when reading dirhb.out! {beta,gamma,energy} not found!"); 
    
    return std::tuple<double,double,double>(beta,gamma,energy);
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
            for( auto const & entry : boost::make_iterator_range( boost::filesystem::directory_iterator("../../output"), {}))
            {
                std::string folder = entry.path().string();

                std::ifstream dirhbout( folder + "/dirhb.out" );

                if( dirhbout.is_open() == false )
                    std::cerr << "Folder: " << folder << ", doesn't contain dirhb.out file!" << std::endl;
                else
                {
                    auto beta_gamma_energy = read_dirhbout( dirhbout );
                    double beta   = std::get<0>(beta_gamma_energy);
                    double gamma  = std::get<1>(beta_gamma_energy);
                    double energy = std::get<2>(beta_gamma_energy);
                    dirhbout.close();

                    // rounding to three digits after decimal point to avoid rounding-error
                    // troubles when insertion and comparison inside a map
                    beta   = std::round( beta   * 1000.0 ) / 1000.0;
                    gamma  = std::round( gamma  * 1000.0 ) / 1000.0;

                    pestable[ std::pair<double,double>(beta,gamma) ] = energy;
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
