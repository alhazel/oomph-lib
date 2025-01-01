// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2025 Matthias Heil and Andrew Hazel
// LIC//
// LIC// This library is free software; you can redistribute it and/or
// LIC// modify it under the terms of the GNU Lesser General Public
// LIC// License as published by the Free Software Foundation; either
// LIC// version 2.1 of the License, or (at your option) any later version.
// LIC//
// LIC// This library is distributed in the hope that it will be useful,
// LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
// LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// LIC// Lesser General Public License for more details.
// LIC//
// LIC// You should have received a copy of the GNU Lesser General Public
// LIC// License along with this library; if not, write to the Free Software
// LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// LIC// 02110-1301  USA.
// LIC//
// LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
// LIC//
// LIC//====================================================================
// Header file for RefineableQElement<1> class
#ifndef OOMPH_REFINEABLE_LINE_ELEMENT_HEADER
#define OOMPH_REFINEABLE_LINE_ELEMENT_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

// oomph-lib headers
#include "binary_tree.h"
#include "macro_element.h"
#include "refineable_elements.h"
#include "Qelements.h"

namespace oomph
{
  // Forward definition for mesh.
  class Mesh;

  //=======================================================================
  /// Refineable version of QElement<1,NNODE_1D>.
  ///
  /// Refinement is performed by binary tree procedures. When the element
  /// is subdivided, the geometry of its sons is established by calls to
  /// their father's \c get_x(...) function which refers to:
  /// - the father element's geometric FE mapping (this is the default)
  ///  or
  /// - a MacroElement's MacroElement::macro_map (if the pointer to the
  /// macro element is non-NULL)
  ///
  /// The class provides a generic RefineableQElement<1>::build() function
  /// which deals with generic isoparametric QElements in which all values
  /// are associated with nodes. The RefineableQElement<1>::further_build()
  /// function provides an interface for any element-specific non-generic
  /// build operations.
  //=======================================================================
  template<>
  class RefineableQElement<1> : public virtual RefineableElement,
                                public virtual LineElementBase
  {
  public:
    /// Shorthand for pointer to an argument-free void member
    /// function of the refineable element
    typedef void (RefineableQElement<1>::*VoidMemberFctPt)();

    /// Constructor: Pass refinement level (default 0 = root)
    RefineableQElement() : RefineableElement()
    {
#ifdef LEAK_CHECK
      LeakCheckNames::RefineableQElement<1> _build += 1;
#endif
    }

    /// Broken copy constructor
    RefineableQElement(const RefineableQElement<1>& dummy) = delete;

    /// Broken assignment operator
    // Commented out broken assignment operator because this can lead to a
    // conflict warning when used in the virtual inheritence hierarchy.
    // Essentially the compiler doesn't realise that two separate
    // implementations of the broken function are the same and so, quite
    // rightly, it shouts.
    /*void operator=(const RefineableQElement<1>&) = delete;*/

    /// Destructor
    virtual ~RefineableQElement()
    {
#ifdef LEAK_CHECK
      LeakCheckNames::RefineableQElement<1> _build -= 1;
#endif
    }

    /// A refineable line element has two sons
    unsigned required_nsons() const
    {
      return 2;
    }

    /// If a neighbouring element has already created a node at a
    /// position corresponding to the local fractional position within the
    /// present element, s_fraction, return a pointer to that node. If
    /// not, return NULL (0). If the node is on a periodic boundary the
    /// flag is_periodic is true, otherwise it will be false.
    Node* node_created_by_neighbour(const Vector<double>& s_fraction,
                                    bool& is_periodic);

    /// If a neighbouring element has already created a node at a
    /// position corresponding to the local fractional position within the
    /// present element, s_fraction, return a pointer to that node. If
    /// not, return NULL (0). If the node is on a periodic boundary the
    /// flag is_periodic is true, otherwise it will be false.
    Node* node_created_by_son_of_neighbour(const Vector<double>& s_fraction,
                                           bool& is_periodic)
    {
      // It is impossible for this situation to arise in meshes
      // containing elements of uniform p-order. This is here so
      // that it can be overloaded for p-refineable elements.
      return 0;
    }

    /// Build the element, i.e. give it nodal positions, apply BCs,
    /// etc. Pointers to any new nodes will be returned in new_node_pt.
    /// If it is open, the positions of the new nodes will be written to
    /// the file stream new_nodes_file.
    virtual void build(Mesh*& mesh_pt,
                       Vector<Node*>& new_node_pt,
                       bool& was_already_built,
                       std::ofstream& new_nodes_file);

    /// Check the integrity of the element: ensure that the position
    /// and values are continuous across the element edges.
    void check_integrity(double& max_error);

    ///  Print corner nodes, using colour
    void output_corners(std::ostream& outfile, const std::string& colour) const;

    /// Pointer to binary tree representation of this element
    BinaryTree* binary_tree_pt()
    {
      return dynamic_cast<BinaryTree*>(Tree_pt);
    }

    /// Pointer to binary tree representation of this element (const version)
    BinaryTree* binary_tree_pt() const
    {
      return dynamic_cast<BinaryTree*>(Tree_pt);
    }

    /// Line elements have no hanging nodes so this is deliberately left empty
    void setup_hanging_nodes(Vector<std::ofstream*>& output_stream) {}

  protected:
    /// Coincidence between son nodal points and father boundaries:
    /// Father_bound[node_1d](jnod_son,son_type) = {L/R/OMEGA}
    static std::map<unsigned, DenseMatrix<int>> Father_bound;

    /// Setup static matrix for coincidence between son nodal points
    /// and father boundaries
    void setup_father_bounds();

    /// Line elements have no hanging nodes so this is deliberately left empty
    void setup_hang_for_value(const int& value_id) {}

    /// Line elements have no hanging nodes so this is deliberately left empty
    void binary_hang_helper(const int& value_id,
                            const int& my_edge,
                            std::ofstream& output_hangfile)
    {
    }
  };


  //=======================================================================
  /// Refineable version of Solid line elements
  //=======================================================================
  template<>
  class RefineableSolidQElement<1> : public virtual RefineableQElement<1>,
                                     public virtual RefineableSolidElement,
                                     public virtual QSolidElementBase
  {
  public:
    /// Constructor, just call the constructor of the RefineableQElement<1>
    RefineableSolidQElement()
      : RefineableQElement<1>(), RefineableSolidElement()
    {
      // Issue a warning about this class
      std::string warning_message =
        "The class RefinableSolidQElement<1> has not been implemented or\n";
      warning_message +=
        "tested. It is safest to assume that all functions do not do what\n";
      warning_message +=
        "they claim to. The `build()' function is deliberately broken.";

      throw OomphLibWarning(
        warning_message,
        "RefineableSolidQElement<1>::RefineableSolidQElement()",
        OOMPH_EXCEPTION_LOCATION);
    }

    /// Broken copy constructor
    RefineableSolidQElement(const RefineableSolidQElement<1>& dummy) = delete;

    /// Broken assignment operator
    /*void operator=(const RefineableSolidQElement<1>&) = delete;*/

    /// Virtual Destructor
    virtual ~RefineableSolidQElement() {}


    /// Final over-ride: Use version in QSolidElementBase
    void set_macro_elem_pt(MacroElement* macro_elem_pt)
    {
      QSolidElementBase::set_macro_elem_pt(macro_elem_pt);
    }

    /// Final over-ride: Use version in QSolidElementBase
    void set_macro_elem_pt(MacroElement* macro_elem_pt,
                           MacroElement* undeformed_macro_elem_pt)
    {
      QSolidElementBase::set_macro_elem_pt(macro_elem_pt,
                                           undeformed_macro_elem_pt);
    }

    /// Use the generic finite difference routine defined in
    /// RefineableSolidElement to calculate the Jacobian matrix
    void get_jacobian(Vector<double>& residuals, DenseMatrix<double>& jacobian)
    {
      RefineableSolidElement::get_jacobian(residuals, jacobian);
    }

    /// Build the element, i.e. give it nodal positions, apply BCs, etc.
    /// Incl. documention into new_nodes_file
    // NOTE: FOR SOME REASON THIS NEEDS TO LIVE IN *.H TO WORK ON INTEL
    void build(Mesh*& mesh_pt,
               Vector<Node*>& new_node_pt,
               bool& was_already_built,
               std::ofstream& new_nodes_file)
    {
      throw OomphLibError("This function has not been implemented yet:",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
  };

} // namespace oomph

#endif
