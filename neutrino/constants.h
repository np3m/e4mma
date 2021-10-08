/// \file Constants.hpp
/// \author jlippuner
/// \since Jul 8, 2014
///
/// \brief
///
///
///
#ifndef EOS_CONSTANTS_HPP_
#define EOS_CONSTANTS_HPP_

#include <cmath>

#ifdef SWIG
%rename(Constants) nuopac::ConstantsStruct;
#endif // SWIG
// we use the workaround http://stackoverflow.com/a/13406244 to make the
// constants show up under the name Constants in Python instead of the constants
// being global constants

namespace nuopac {

  ///
  /// Structure containing a number of useful physical constants
  ///
  struct ConstantsStruct {
    /// GFermi/(hbar c)^3 in [MeV^-2] 
    static constexpr double GfMeV = 1.1663787e-11;

    /// Sin(\theta_C) 
    static constexpr double sinThetaC = 0.231;
  
    /// DWeak and Dweak 
    static constexpr double DWeak = 0.756;
    static constexpr double FWeak = 0.477;

    /// Sin^2(\theta_W) 
    static constexpr double sin2ThetaW = 0.21320; 
  
    /// pi
    static constexpr double Pi = 3.1415926535897932385;
  
    /// Hbar * c in [MeV fm]
    static constexpr double HBCFmMeV = 197.3269788; 
   
    /// erg per MeV http://physics.nist.gov/cgi-bin/cuu/Value?tevj
    static constexpr double ErgPerMeV = 1.602176565E-6;

    /// c http://physics.nist.gov/cgi-bin/cuu/Value?c
    static constexpr double SpeedOfLightInCmPerSec = 2.99792458E10;

    /// hbar http://physics.nist.gov/cgi-bin/cuu/Value?hbarev
    static constexpr double ReducedPlanckConstantInMeVSec = 6.58211928E-22;

    /// k_B http://physics.nist.gov/cgi-bin/cuu/Value?tkev
    static constexpr double BoltzmannConstantInMeVPerGK = 8.6173324E-2;

    static constexpr double BoltzmannConstantsInErgPerK =
      BoltzmannConstantInMeVPerGK * 1.0E-9 * ErgPerMeV;

    /// N_A http://physics.nist.gov/cgi-bin/cuu/Value?na
    static constexpr double AvogadroConstantInPerGram = 6.02214129E23;

    /// 2 pi hbar^2 c^2
    static constexpr double TwoPiHbar2C2 = (2.0 * Pi
                                            * ReducedPlanckConstantInMeVSec * ReducedPlanckConstantInMeVSec
                                            * SpeedOfLightInCmPerSec * SpeedOfLightInCmPerSec);

    //  /// 2 pi hbar^2 c^2 N_A^{2/3}
    //#if defined(__ICC) || defined(__INTEL_COMPILER)
    //  // we cannot use pow with the Intel compiler, it's a bug:
    //  // https://software.intel.com/en-us/forums/topic/484936
    //  static constexpr double TwoPiHbar2C2NA23 = TwoPiHbar2C2 * 7.13127680E15;
    //
    //#else
    //  static constexpr double TwoPiHbar2C2NA23 = TwoPiHbar2C2
    //      * pow(AvogadroConstantInPerGram, 2.0 / 3.0);
    //#endif

    /// k_B / (2 pi hbar^2 c^2)
    static constexpr double BoltzmannConstantDivBy2PiHbar2C2InPerMeVPerGKPerCm2 =
      BoltzmannConstantInMeVPerGK / TwoPiHbar2C2;

    /// k_B / (2 pi hbar^2 c^2 N_A^{2/3})
    //static constexpr double InverseRateFactor = BoltzmannConstantInMeVPerGK
    //    / TwoPiHbar2C2NA23;

    /// m_p http://physics.nist.gov/cgi-bin/cuu/Value?mpc2mev
    static constexpr double ProtonMassInMeV = 938.272046;
  
    static constexpr double ProtonMassInFm = ProtonMassInMeV / HBCFmMeV;

    /// m_e http://physics.nist.gov/cgi-bin/cuu/Value?mec2mev
    static constexpr double ElectronMassInMeV = 0.510998928;
  
    static constexpr double ElectronMassInFm = 0.510998928 / HBCFmMeV;

    /// m_p + m_e
    static constexpr double ProtonMassPlusElectronMassInMeV =
      ProtonMassInMeV + ElectronMassInMeV;

    /// m_n http://physics.nist.gov/cgi-bin/cuu/Value?mnc2mev
    static constexpr double NeutronMassInMeV = 939.565379;
  
    static constexpr double NeutronMassInFm = 939.565379 / HBCFmMeV;
 
    /// Elementary charge (unitless when hbar = c = 1)
    static constexpr double ElementaryChargeSquared = 1.4299764/HBCFmMeV;
   
  };

  typedef ConstantsStruct Constants;

}
#endif // EOS_CONSTANTS_HPP_
