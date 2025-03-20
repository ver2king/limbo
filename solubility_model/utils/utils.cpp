#include "utils.hpp"

namespace CONST
{   
    double const minDivTol = 1E-12;
    double const gasR = 83.1447;
    double const refP = 1.; // bar
    double const saltV = 2.; 

    double const floatInf = std::numeric_limits<float>::inifinity();
    double const doubleInf = std::numeric_limits<double>::inifinity();
}

namespace UNIT
{
    int ConvertPressureUnit(std::pair<std::string, std::string> unitPair)
    {
        if ( unitPair.first.compare("MPa") == 0 && unitPair.second.compare("Bar") == 0 )
        { return 0; }
        else if ( unitPair.first.compare("KPa") == 0 && unitPair.second.compare("Bar") == 0  )
        { return 1; }
        else if ( unitPair.first.compare("GPa") == 0 && unitPair.second.compare("Bar") == 0  )
        { return 2; }
        else if ( unitPair.first.compare("Pa") == 0 && unitPair.second.compare("Bar") == 0  )
        { return 3; }
        else if ( unitPair.first.compare("Psi") == 0 && unitPair.second.compare("Bar") == 0 )
        { return 4; }
        else { return -1; }; // No unit conversion
        return -1;
    };

    int ConvertTemperatureUnit(std::pair<std::string, std::string> unitPair)
    {
        if ( unitPair.first.compare("C") == 0 && unitPair.second.compare("K") == 0 )
        { return 0; }
        else if ( unitPair.first.compare("F") == 0 && unitPair.second.compare("K") == 0  )
        { return 1; }
        else if ( unitPair.first.compare("F") == 0 && unitPair.second.compare("C") == 0  )
        { return 2; }
        else if ( unitPair.first.compare("K") == 0 && unitPair.second.compare("C") == 0  )
        { return 3; }
        else { return -1; }; // No unit conversion
        return -1;
    };

    double ConvertPressure(double P, std::pair<std::string, std::string> unitPair)
    {
        double defP = std::numeric_limits<double>::quiet_NaN();
        int convertCode = ConvertPressureUnit(unitPair);
        switch (convertCode)
        {
        case 0:
            return 10. * P;
        case 1:
            return 1e2 * P;
        case 2:
            return 10000. * P;
        case 3:
            return 1e-5 * P;
        case 4:
            return 0.0689475729* P;
        case -1:
            return P;
        default:
            return defP;
        };
    };

    double ConvertTemperature(double T, std::pair<std::string, std::string> unitPair)
    {
        double defT = std::numeric_limits<double>::quiet_NaN();
        int convertCode = ConvertTemperatureUnit(unitPair);
        switch (convertCode)
        {
        case 0:
            return 273.15 + T;
        case 1:
            return ( T - 32.) / 1.8 + 273.15;
        case 2:
            return ( T - 32.) / 1.8;
        case 3:
            return T - 273.15;
        case -1:
            return T;
        default:
            return defT;
        };
    };
}

namespace UTILS
{

    Tensor2DFloat64 ReadMatrixFile(std::ifstream &f, int &row, int &col)
    {   

        Tensor2DFloat64 val(row, std::vector<double>(col, 0.));
        Tensor1DFloat64 val_;
        String data;
        char delim = ' ';
        //
        if (f.is_open())
        {
        while ( std::getline(f, data) )
        {   
            String d;
            std::stringstream line_data(data);
            while( std::getline(line_data, d, delim) )
            {   
                val_.push_back( std::stod(d) );
            };
        };
        f.close();
        } else
        {
        printf("File not opened.");
        };
        //
        for ( Index ir = 0; ir < row; ++ ir )
        {
            for ( Index ic = 0; ic < col; ++ ic )
            {   
                Index i = ir * col + ic;
                val[ir][ic] = val_[i];
            };
        };
        return val;
    };

    Tensor2DFloat64 ReadMatrixFileMultiDelimiters(std::ifstream &f, int &row, int &col)
    {       
    std::regex re("\t");
    Tensor2DFloat64 val(row, std::vector< double > ( col, 0. ));
    Tensor1DFloat64 val_;
    String data;
    //
    if (f.is_open())
    {
    int val_Idx = 0;
    while ( std::getline(f, data) )
    {   
        std::sregex_token_iterator first{data.begin(), data.end(), re, -1}; 
        std::sregex_token_iterator last; 
        std::vector<std::string> tokens{first, last};
        for ( std::string str : tokens )
        { 
            if ( str.compare("\t") != 0 ) { val_.push_back( std::stod( str ) ); };
        }
    };
    f.close();
    } else
    {
    printf("File not opened.");
    };
    //
    for ( int ir = 0; ir < row; ++ ir )
    {
        for ( int ic = 0; ic < col; ++ ic )
        {   
            int i = ir * col + ic;
            val[ir][ic] = val_[i];
        };
    };
    return val;
    };

    double LinearScaler( double inputVal, double lowerBound, double higherBound )
    {
        /// Scale an input value in [0, 1] linearly to [lowerBound, higherBound]
        double scaledVal = lowerBound + inputVal * ( higherBound - lowerBound );
        return scaledVal;
    };

    Tensor1DFloat64 CubicEquationSolver(Tensor1DFloat64 cubicCoeffs)
    {   
        //
        double a=cubicCoeffs[0];
        double b=cubicCoeffs[1];
        double c=cubicCoeffs[2];
        double d=cubicCoeffs[3];
        b /= a;
        c /= a;
        d /= a;
        //
        double disc, q, r, dum1, s, t, term1, r13;
        q = (3.0*c - (b*b))/9.0;
        r = -(27.0*d) + b*(9.0*c - 2.0*(b*b));
        r /= 54.0;
        disc = q*q*q + r*r;
        term1 = b/3.0;
        if (disc > 0) { 
        // one root real, two are complex
        s = r + sqrt(disc);
        s = ((s < 0) ? -pow(-s, (1.0/3.0)) : pow(s, (1.0/3.0)));
        t = r - sqrt(disc);
        t = ((t < 0) ? -pow(-t, (1.0/3.0)) : pow(t, (1.0/3.0)));
        double x1=-term1 + s + t;
        Tensor1DFloat64 roots(3, 0.);
        roots[0]=x1; 
        roots[1]=std::numeric_limits<double>::infinity();
        roots[2]=std::numeric_limits<double>::infinity();
        return roots; 
        } 
        else if (disc == 0) { 
        // All roots real, at least two are equal.
        r13 = ((r < 0) ? -pow(-r,(1.0/3.0)) : pow(r,(1.0/3.0)));
        double x1 = -term1 + 2.0*r13;
        double x2 = -r13 + term1;
        Tensor1DFloat64 roots(3, 0.);
        roots[0]=x1; roots[1]=x2; roots[2]=x2; 
        return roots; 
        } else { 
        // Only option left is that all roots are real and unequal (to get here, q < 0)
        q = -q;
        dum1 = q*q*q;
        dum1 = acos(r/sqrt(dum1));
        r13 = 2.0*sqrt(q);
        double x1 = -term1 + r13*cos(dum1/3.0);
        double x2 = -term1 + r13*cos((dum1 + 2.0*M_PI)/3.0);
        double x3 = -term1 + r13*cos((dum1 + 4.0*M_PI)/3.0);
        Tensor1DFloat64 roots(3, 0.);
        roots[0]=x1; roots[1]=x2; roots[2]=x3; 
        return roots;
        };
    };
}

