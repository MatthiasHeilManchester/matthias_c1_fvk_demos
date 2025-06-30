#ifndef OOMPH_NEW_INTEGRATION_HEADER
#define OOMPH_NEW_INTEGRATION_HEADER

// Include the Integral class
#include "generic/integral.h"

namespace oomph
{
  //=========================================================
  /// Class for integration rule for integrands weighted by
  /// x^(-1/2) singularity
  ///
  /// Empty -- just establishes the template parameters.
  //=========================================================
  template<unsigned DIM, unsigned NPTS_1D>
  class SingularIntegration
  {
  };


  //=========================================================
  /// 1D x^(-1/2) weighted integrand integration class.
  /// One integration point. This integration scheme can
  /// integrate up to first-order polynomials multiplied by
  /// x^(-1/2) exactly.
  //=========================================================
  template<>
  class SingularIntegration<1, 1> : public Integral
  {
  private:
    /// Number of integration points in the scheme
    static const unsigned Npts = 1;
    /// Array to hold weight and knot points (defined in cc file)
    static const double Knot[1][1], Weight[1];

  public:
    /// Default constructor (empty)
    SingularIntegration(){};

    /// Broken copy constructor
    SingularIntegration(const SingularIntegration& dummy) = delete;

    /// Broken assignment operator
    void operator=(const SingularIntegration&) = delete;

    /// Number of integration points of the scheme
    unsigned nweight() const
    {
      return Npts;
    }

    /// Return coordinate x[j] (j=0) of integration point i
    double knot(const unsigned& i, const unsigned& j) const
    {
      return Knot[i][j];
    }

    /// Return weight of integration point i
    double weight(const unsigned& i) const
    {
      return Weight[i];
    }
  };


  //=========================================================
  /// 1D x^(-1/2) weighted integrand integration class.
  /// Two integration points. This integration scheme can
  /// integrate up to third-order polynomials weighted by
  /// x^(-1/2) exactly and is therefore a suitable "full"
  /// integration scheme for linear (two-node) elements in
  /// which the highest-order polynomial is quadratic.
  //=========================================================
  template<>
  class SingularIntegration<1, 2> : public Integral
  {
  private:
    /// Number of integration points in the scheme
    static const unsigned Npts = 2;
    /// Array to hold weight and knot points (defined in cc file)
    static const double Knot[2][1], Weight[2];

  public:
    /// Default constructor (empty)
    SingularIntegration(){};

    /// Broken copy constructor
    SingularIntegration(const SingularIntegration& dummy) = delete;

    /// Broken assignment operator
    void operator=(const SingularIntegration&) = delete;

    /// Number of integration points of the scheme
    unsigned nweight() const
    {
      return Npts;
    }

    /// Return coordinate x[j] (j=0) of integration point i
    double knot(const unsigned& i, const unsigned& j) const
    {
      return Knot[i][j];
    }

    /// Return weight of integration point i
    double weight(const unsigned& i) const
    {
      return Weight[i];
    }
  };


  //=========================================================
  /// 1D x^(-1/2) weighted integrand integration class.
  /// Three integration points. This integration scheme can
  /// integrate up to fifth-order polynomials weighted by
  /// x^(-1/2) exactly and is therefore a suitable "full"
  /// integration scheme for quadratic (three-node) elements in
  /// which the highest-order polynomial is fourth order.
  //=========================================================
  template<>
  class SingularIntegration<1, 3> : public Integral
  {
  private:
    /// Number of integration points in the scheme
    static const unsigned Npts = 3;
    /// Array to hold weight and knot points (defined in cc file)
    static const double Knot[3][1], Weight[3];

  public:
    /// Default constructor (empty)
    SingularIntegration(){};

    /// Broken copy constructor
    SingularIntegration(const SingularIntegration& dummy) = delete;

    /// Broken assignment operator
    void operator=(const SingularIntegration&) = delete;

    /// Number of integration points of the scheme
    unsigned nweight() const
    {
      return Npts;
    }

    /// Return coordinate x[j] (j=0) of integration point i
    double knot(const unsigned& i, const unsigned& j) const
    {
      return Knot[i][j];
    }

    /// Return weight of integration point i
    double weight(const unsigned& i) const
    {
      return Weight[i];
    }
  };


  //=========================================================
  /// 1D x^(-1/2) weighted integrand integration class.
  /// Four integration points. This integration scheme can
  /// integrate up to seventh-order polynomials weighted by
  /// x^(-1/2) exactly and is therefore a suitable "full"
  /// integration scheme for cubic (four-node) elements in
  /// which the highest-order polynomial is sixth order.
  //=========================================================
  template<>
  class SingularIntegration<1, 4> : public Integral
  {
  private:
    /// Number of integration points in the scheme
    static const unsigned Npts = 4;
    /// Array to hold weight and knot points (defined in cc file)
    static const double Knot[4][1], Weight[4];

  public:
    /// Default constructor (empty)
    SingularIntegration(){};

    /// Broken copy constructor
    SingularIntegration(const SingularIntegration& dummy) = delete;

    /// Broken assignment operator
    void operator=(const SingularIntegration&) = delete;

    /// Number of integration points of the scheme
    unsigned nweight() const
    {
      return Npts;
    }

    /// Return coordinate x[j] (j=0) of integration point i
    double knot(const unsigned& i, const unsigned& j) const
    {
      return Knot[i][j];
    }

    /// Return weight of integration point i
    double weight(const unsigned& i) const
    {
      return Weight[i];
    }
  };

  //=========================================================
  /// 1D x^(-1/2) weighted integrand integration class.
  /// Five integration points. This integration scheme can
  /// integrate up to ninth-order polynomials weighted by
  /// x^(-1/2) exactly and is therefore a suitable "full"
  /// integration scheme for quartic (five-node) elements in
  /// which the highest-order polynomial is eighth order.
  //=========================================================
  template<>
  class SingularIntegration<1, 5> : public Integral
  {
  private:
    /// Number of integration points in the scheme
    static const unsigned Npts = 5;
    /// Array to hold weight and knot points (defined in cc file)
    static const double Knot[5][1], Weight[5];

  public:
    /// Default constructor (empty)
    SingularIntegration(){};

    /// Broken copy constructor
    SingularIntegration(const SingularIntegration& dummy) = delete;

    /// Broken assignment operator
    void operator=(const SingularIntegration&) = delete;

    /// Number of integration points of the scheme
    unsigned nweight() const
    {
      return Npts;
    }

    /// Return coordinate x[j] (j=0) of integration point i
    double knot(const unsigned& i, const unsigned& j) const
    {
      return Knot[i][j];
    }

    /// Return weight of integration point i
    double weight(const unsigned& i) const
    {
      return Weight[i];
    }
  };


  //=========================================================
  /// 1D x^(-1/2) weighted integrand integration class.
  /// Six integration points. This integration scheme can
  /// integrate up to eleventh-order polynomials weighted by
  /// x^(-1/2) exactly and is therefore a suitable "full"
  /// integration scheme for quintic (six-node) elements in
  /// which the highest-order polynomial is tenth order.
  //=========================================================
  template<>
  class SingularIntegration<1, 6> : public Integral
  {
  private:
    /// Number of integration points in the scheme
    static const unsigned Npts = 6;
    /// Array to hold weight and knot points (defined in cc file)
    static const double Knot[6][1], Weight[6];

  public:
    /// Default constructor (empty)
    SingularIntegration(){};

    /// Broken copy constructor
    SingularIntegration(const SingularIntegration& dummy) = delete;

    /// Broken assignment operator
    void operator=(const SingularIntegration&) = delete;

    /// Number of integration points of the scheme
    unsigned nweight() const
    {
      return Npts;
    }

    /// Return coordinate x[j] (j=0) of integration point i
    double knot(const unsigned& i, const unsigned& j) const
    {
      return Knot[i][j];
    }

    /// Return weight of integration point i
    double weight(const unsigned& i) const
    {
      return Weight[i];
    }
  };
} // namespace oomph

#endif