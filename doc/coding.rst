
Coding conventions
==================

General guidelines
++++++++++++++++++

In general all code in the project should follow the following
conventions.  Sometimes a few exceptions may be acceptable if it is
not possible/practical to follow these rules.

#. Use only 80 columns. Break lines if needed.

#. Use lowercase first and (perhaps) underscores for function names:
   copy_to()

#. Classes are capitalized: Element, MyClass

#. Use lowercase first and no underscores for 'regular' public
   attributes: vt0, isNonlinear

#. For attributes added to arbitrary instances use a unique prefix
   that ends with and underscore to avoid conflict with internal
   attributes. Example: ``nD_`` (used in nodal.py)

#. Use underscore first for private attributes: _privateVar 

#. module names: first letter lowercase and no underscores: circuit,
   paramset


Documentation
+++++++++++++

Keep doctrings for modules up to date. Consider writing some
preliminary documentation before implementing a new feature.  That
documentation will be included in the module reference. The User's
Manual entry for device models is automatically generated from the doc
string of the ``Device`` class. Similarly, for analyses use the
doctring from the ``Analysis`` class in the corresponding module.

Version Numbers
+++++++++++++++

Starting with version 0.4.1, ``.dev`` is appended at the end of a
version number while the version evolves, *i.e.*, changes are commited
to that version. Once a version is deemed to be released, the ``.dev``
part should be stripped and the commit tagged with
``version-0.x.y``. An increment in ``x`` indicates that a major
feature has been added (such as a new analysis type, major re-write of
code, etc.). Otherwise, ``y`` is incremented.

Version 1.0 will be released when the following features are
implemented:

  * Basic analyses are efficient: DC, AC, Transient
 
  * Frequency-defined devices are supported in all analyses

  * Time-delay interface is stable and implemented for all analyses

  * Netlist syntax and interfaces in general are mostly stable
