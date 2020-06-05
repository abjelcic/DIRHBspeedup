#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <exception>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>


enum class EDF { DDPC1, DDME2 };
std::map< EDF , std::string > EDFnames = { { EDF::DDPC1 , "DD-PC1" },
                                           { EDF::DDME2 , "DD-ME2" }
                                          };

using REAL = double;
int generate_dirhbdat( std::string nucleus,
                       std::size_t A,

                       std::size_t n0f,
                       std::size_t n0b,

                       EDF parameterset,

                       REAL beta0, REAL gamma0,
                       REAL betai, REAL gammai,
                       REAL betac, REAL gammac,

                       REAL cquad,

                       std::string  path      = "./",

                       bool         icstr     = true,
                       bool         inin_wel  = true,
                       bool         inin_del  = true,
                       REAL         gap_n     = REAL(1),
                       REAL         gap_p     = REAL(1)
                      )
{

	std::ofstream dirhbdat( path+"dirhb.dat" , std::ios::out );

	if( dirhbdat.is_open() == false )
		return -1;

	try
	{
		dirhbdat.precision(3);

		dirhbdat << "n0f,n0b  ="
				 << std::setw(4) << n0f
				 << " "
				 << std::setw(4) << n0b
				 << "                ! number of oscillator shells(F,B)"
				 << std::endl;

		dirhbdat << "beta0    ="
				 << std::setw(9) << std::fixed << std::showpos << beta0 << std::noshowpos
				 << "                ! beta-deformation parameter of the basis"
				 << std::endl;

		dirhbdat << "gamma0   ="
				 << std::setw(9) << std::fixed << gamma0
				 << "                ! gamma-deformation parameter of the basis"
				 << std::endl;

		dirhbdat << "betai    ="
				 << std::setw(9) << std::fixed << std::showpos << betai << std::noshowpos
				 << "                ! beta-deformation  of the start potential"
				 << std::endl;

		dirhbdat << "gammai   ="
				 << std::setw(9) << std::fixed << gammai
				 << "                ! gamma-deformation of the start potential"
				 << std::endl;

		dirhbdat << "inin     ="
				 << std::setw(5) << ( inin_wel ? 1 : 0 )
				 << std::setw(5) << ( inin_del ? 1 : 0 )
				 << "               ! initialization"
				 << std::endl;

		dirhbdat.precision(2);
		dirhbdat << "init. gap="
				 << std::setw(8) << std::fixed << gap_n
				 << std::setw(8) << std::fixed << gap_p
				 << std::endl;

		dirhbdat << nucleus.substr(0,2) << " " << std::setw(3) << A
				 << "                             ! nucleus under consideration"
				 << std::endl;

		dirhbdat << "c-------------------------------------------------------------------" << std::endl;


		dirhbdat << "Force    =  "
				 << EDFnames[parameterset]
				 << "                 ! Parameterset of the Lagrangian"
				 << std::endl;

		dirhbdat << "c-------------------------------------------------------------------" << std::endl;

		dirhbdat << "icstr    ="
				 << std::setw(5) << ( icstr ? 1 : 0 )
				 << "                    ! constraint none (0) quadratic (1)"
				 << std::endl;

		dirhbdat.precision(3);
		dirhbdat << "betac    ="
				 << std::setw(9) << std::fixed << std::showpos << betac << std::noshowpos
				 << "                ! constrained value of beta"
				 << std::endl;

		dirhbdat << "gammac   ="
				 << std::setw(9) << std::fixed << gammac
				 << "                ! constrained value of gamma"
				 << std::endl;

		dirhbdat << "cqad     ="
				 << std::setw(9) << std::fixed << cquad
				 << "                ! Spring constant for constraint"
				 << std::endl;

		dirhbdat << "c-------------------------------------------------------------------" << std::endl;

	}
	catch( ... )
	{
		std::cerr << "Something went terribly wrong..." << std::endl;
		return -2;
	}




	dirhbdat.close();
	return 0;
}















int main(int argc, char *argv[])
{

    std::ifstream pesdat( "../../pes.dat" );
	if( pesdat.is_open() == false )
    {
        std::cerr << "pes.dat opening failed!" << std::endl;
        return -1;
    }


    std::string     nucleus;
    std::size_t     A;
    std::size_t     n0f;
    std::size_t     n0b;
    EDF             parameterset;
    double          beta_start;
    double          beta_end;
    double          delta_beta;
    double          gamma_start;
    double          gamma_end;
    double          delta_gamma;
    double          cquad;

    try
    {
        std::string line;
        while( std::getline(pesdat,line) )
        {
            boost::trim_right(line);
            boost::trim_left (line);

            if( line == "#Element name" )
            {
                std::getline(pesdat,line);
                std::stringstream(line) >> nucleus;
                continue;
            }
            if( line == "#Number of nucleons A" )
            {
                std::getline(pesdat,line);
                std::stringstream(line) >> A;
                continue;
            }
            if( line == "#Parameterset of the Lagrangian (DD-PC1 or DD-ME2)" )
            {
                std::getline(pesdat,line);
                std::string parname;
                std::stringstream(line) >> parname;

                bool found = false;
                for(auto kv : EDFnames)
                    if( kv.second == parname )
                    {
                        parameterset = kv.first;
                        found = true;
                        break;
                    }
                if( found == false )
                {
                    std::cerr << "Parameterset of the Lagrangian invalid!" << std::endl;
                    return -2;
                }

                continue;
            }
            if( line == "#Number of shells (n0f)" )
            {
                std::getline(pesdat,line);
                std::stringstream(line) >> n0f;
                continue;
            }
            if( line == "#Number of shells for bosons (n0b)" )
            {
                std::getline(pesdat,line);
                std::stringstream(line) >> n0b;
                continue;
            }
            if( line == "#Spring constant for constraint (cqad)" )
            {
                std::getline(pesdat,line);
                std::stringstream(line) >> cquad;
                continue;
            }
            if( line == "#Hill-Wheeler beta sweep (beta_start, beta_end, delta_beta)" )
            {
                std::getline(pesdat,line);
                std::stringstream(line) >> beta_start;
                std::getline(pesdat,line);
                std::stringstream(line) >> beta_end;
                std::getline(pesdat,line);
                std::stringstream(line) >> delta_beta;
                continue;
            }
            if( line == "#Hill-Wheeler gamma sweep (gamma_start, gamma_end, delta_gamma)" )
            {
                std::getline(pesdat,line);
                std::stringstream(line) >> gamma_start;
                std::getline(pesdat,line);
                std::stringstream(line) >> gamma_end;
                std::getline(pesdat,line);
                std::stringstream(line) >> delta_gamma;
                continue;
            }

        }
    }
    catch( ... )
    {
		std::cerr << "Something went terribly wrong..." << std::endl;
		return -3;
    }

    pesdat.close();

/**
    std::cout << nucleus                << std::endl;
    std::cout << A                      << std::endl;
    std::cout << n0f                    << std::endl;
    std::cout << n0b                    << std::endl;
    std::cout << EDFnames[parameterset] << std::endl;
    std::cout << beta_start             << std::endl;
    std::cout << beta_end               << std::endl;
    std::cout << delta_beta             << std::endl;
    std::cout << gamma_start            << std::endl;
    std::cout << gamma_end              << std::endl;
    std::cout << delta_gamma            << std::endl;
    std::cout << cquad                  << std::endl;
**/






    try
    {
        if( !boost::filesystem::exists( "../../output" ) )
            boost::filesystem::create_directory("../../output");

        for( auto beta = beta_start; beta <= beta_end+1.e-8; beta += delta_beta )
            for( auto gamma = gamma_start; gamma <= gamma_end+1.e-8; gamma += delta_gamma )
            {

                std::ostringstream oss;
                oss.precision(3);
                oss << std::fixed << std::showpos;
                oss << "b" << beta << "g" << gamma;

                std::string dirname = "../../output/" + nucleus + std::to_string(A) + "_" + oss.str();

                if( !boost::filesystem::exists(dirname) )
                {
                    boost::filesystem::create_directory(dirname);
                    boost::filesystem::copy_file( "../dirhbt/run" , dirname + "/run" + oss.str() );

                    auto report = generate_dirhbdat(nucleus, A, n0f, n0b, parameterset,
                                                    beta, gamma, beta, gamma, beta, gamma, cquad,
                                                    dirname+"/" );
                    if( report != 0 )
                        throw std::runtime_error("generate_dirhbdat failed with code " + std::to_string(report));

                }

            }
    }
    catch( std::exception const & e )
    {
        std::cout << "Exception: " << e.what() << std::endl;
    }




	return 0;
}


