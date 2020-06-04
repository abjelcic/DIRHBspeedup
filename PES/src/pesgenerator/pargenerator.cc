#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <boost/algorithm/string.hpp>

int generate_dirhbpar( std::size_t n0f,
                       std::size_t n0b,
                       std::size_t NGH,
                       std::string path = "./"
                      )
{

	std::ofstream dirhbpar( path+"dirhb.par" );

	if( dirhbpar.is_open() == false )
		return -1;

	try
	{
        dirhbpar << "c-----------------------------------------------------------------------    " << std::endl
                 << "c     Parameter file for:                                                   " << std::endl
                 << "c     n0f =" << std::setw(4) << n0f << "  n0b =" << std::setw(4) << n0b << "" << std::endl
                 << "c-----------------------------------------------------------------------    " << std::endl
                 << "c                                                                           " << std::endl
                 << "c---- maximal number for GFV                                                " << std::endl
                 << "      parameter (   igfv =       100 )                                      " << std::endl
                 << "c                                                                           " << std::endl
                 << "c---- number of gauss-meshpoints                                            " << std::endl
                 << "      parameter (    ngh =" << std::setw(10) << NGH << " )                  " << std::endl
                 << "                                                                            " << std::endl
                 << "c---- maximal number of (k,parity)-blocks                                   " << std::endl
                 << "      parameter (    nbx =        2 )                                       " << std::endl
                 << "c                                                                           " << std::endl
                 << "c---- maximal oscillator quantum number for large components                " << std::endl
                 << "      parameter (   n0fx =" << std::setw(3) << n0f << " )                   " << std::endl
                 << "      parameter (   kmaxp= n0fx/2  )                                        " << std::endl
                 << "      parameter (   n0fp = (kmaxp+1)*(kmaxp+2)*(4*kmaxp+3)/6)               " << std::endl
                 << "                                                                            " << std::endl
                 << "      parameter (   kmaxn= n0fx/2 -1+mod(n0fx,2) )                          " << std::endl
                 << "      parameter (   n0fn = (kmaxn+1)*(kmaxn+2)*(4*kmaxn+9)/6)               " << std::endl
                 << "      parameter (   nfx = max(n0fp,n0fn)  )                                 " << std::endl
                 << "                                                                            " << std::endl
                 << "c---- maximal oscillator quantum number for small componente                " << std::endl
                 << "      parameter (   n0gx =       n0fx+1 )                                   " << std::endl
                 << "      parameter (   kmaxp1= n0gx/2  )                                       " << std::endl
                 << "      parameter (   n0gp = (kmaxp1+1)*(kmaxp1+2)*(4*kmaxp1+3)/6)            " << std::endl
                 << "                                                                            " << std::endl
                 << "      parameter (   kmaxn1= n0gx/2 -1+mod(n0gx,2) )                         " << std::endl
                 << "      parameter (   n0gn = (kmaxn1+1)*(kmaxn1+2)*(4*kmaxn1+9)/6)            " << std::endl
                 << "      parameter (   ngx = max(n0gp,n0gn)  )                                 " << std::endl
                 << "                                                                            " << std::endl
                 << "      parameter (   ndx = max(nfx,ngx)  )                                   " << std::endl
                 << "                                                                            " << std::endl
                 << "c---- max. number of eigenstates for protons or neutrons                    " << std::endl
                 << "      parameter (    nkx =       n0gp+n0gn )                                " << std::endl
                 << "c                                                                           " << std::endl
                 << "c---- max. number of all levels for protons or neutrons                     " << std::endl
                 << "      parameter (    ntx =      n0fp+n0fn+n0gp+n0gn )                       " << std::endl
                 << "                                                                            " << std::endl
                 << "c---- maximal oscillator quantum number for bosons                          " << std::endl
                 << "      parameter (   n0bx =" << std::setw(10) << n0b << " )                  " << std::endl
                 << "c---- maximal dimension for boson-matrix                                    " << std::endl
                 << "      parameter (   nobx =  (n0bx+2)*(n0bx+4)*(n0bx+6)/48      )            " << std::endl
                 << "      parameter (   nox1 = (2*n0gx+2)*(2*n0gx+4)*(2*n0gx+6)/48 )            " << std::endl
                 << "c                                                                           " << std::endl
                 << "c                                                                           " << std::endl
                 << "c---- maximal n for Hermite polynomials                                     " << std::endl
                 << "      parameter (   nmax =" << std::setw(10) << std::max(n0f+1,n0b) << " )  " << std::endl
                 << "c                                                                           " << std::endl
                 << "c-----------------------------------------------------------------------    " << std::endl
                 << "c                                                                           " << std::endl
                 << "      parameter (nb2x  = nbx+nbx)                                           " << std::endl
                 << "      parameter (nhx   = nfx+ngx)                                           " << std::endl
                 << "      parameter (nddx  = ndx*ndx)                                           " << std::endl
                 << "      parameter (nhfbx = nhx+nhx)                                           " << std::endl
                 << "      parameter (nhfbqx= nhfbx*nhfbx)                                       " << std::endl
                 << "      parameter (nfgx  = nfx*ngx)                                           " << std::endl
                 << "      parameter (nhhx  = nhx*nhx)                                           " << std::endl
                 << "      parameter (noox  = nobx*nobx)                                         " << std::endl
                 << "      parameter (ngauss= (ngh+1)*(ngh+1)*(ngh+1))                           " << std::endl
                 << "      parameter (kmax  = n0fx+1)                                            " << std::endl
                 << "      parameter (kmax2 = 2*kmax)                                            " << std::endl
                 << "                                                                            " << std::endl
                 << "      parameter (lwork=1+6*nhfbx+2*nhfbx**2)                                " << std::endl
                 << "      parameter (liwork=3+5*nhfbx)                                          " << std::endl
                 << "c-----Pairing interaction                                                   " << std::endl
                 << "      parameter (npmax = 2*n0fx)                                            " << std::endl
                 << "      parameter (np2   = npmax/2+1)                                         " << std::endl
                 << "      parameter (nsepx = (npmax+2)*(npmax+4)*(npmax+6)/48)                  " << std::endl
                 << "                                                                            " << std::endl
                 << "      parameter(mvfp = n0fp*(n0fp+1)/2 )                                    " << std::endl
                 << "      parameter(mvfn = n0fn*(n0fn+1)/2 )                                    " << std::endl
                 << "      parameter(mvfx = mvfp+mvfn )                                          " << std::endl
                 << "                                                                            " << std::endl
                 << "      parameter(mvgp = n0gp*(n0gp+1)/2)                                     " << std::endl
                 << "      parameter(mvgn = n0gn*(n0gn+1)/2)                                     " << std::endl
                 << "      parameter(mvgx = mvgp+mvgn)                                           " << std::endl
                 << "                                                                            " << std::endl
                 << "      parameter(mvtp= (n0fp+n0gn)*(n0fp+n0gn+1)/2 )                         " << std::endl
                 << "      parameter(mvtn= (n0fn+n0gp)*(n0fn+n0gp+1)/2 )                         " << std::endl
                 << "      parameter(mvtx= mvtp+mvtn)                                            " << std::endl
                 << "c-----Coulomb interaction                                                   " << std::endl
                 << "c                                                                           " << std::endl
                 << "      parameter (ndcoul = 160)                                              " << std::endl
                 << "      parameter (ndneta = 159)                                              " << std::endl
                 << "      parameter (furmax = 0.25)                                             " << std::endl
                 << "      parameter (boucou = 40.0)                                             " << std::endl
                 << "      parameter (ndkart = 3)                                                " << std::endl
                 << "      parameter (nupols = 6)                                                " << std::endl;


	}
	catch( ... )
	{
		std::cerr << "Something went terribly wrong..." << std::endl;
		return -2;
	}




	dirhbpar.close();
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



    std::size_t n0f;
    std::size_t n0b;

    std::string line;
    while( std::getline(pesdat,line) )
    {
        boost::trim_right(line);
        boost::trim_left (line);

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
    }

    pesdat.close();


    std::size_t NGH = std::max( 48ul , std::max( 2*(n0f+1) , 2*n0b ) );
    if( generate_dirhbpar(n0f, n0b, NGH, "../dirhbt/") )
    {
        std::cerr << "generate_dirhbpar failed!" << std::endl;
        return -2;
    }


	return 0;
}


