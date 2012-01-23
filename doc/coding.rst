
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

Keep doctrings for modules up to date. That documentation will be
included in the module reference. The user manual entry for device
models is automatically generated from the doc string of the
``Device`` class. Similarly, for analyses use the doctring from the
``Analysis`` class in the corresponding module.

Version Numbers
+++++++++++++++

The initial release of the program was version 0.1.1 and had no basic
analysis implemented. Starting from version 0.2.dev, versions are finished
with a ``.dev`` tag to indicate that the program is still being
developed and interfaces may change.  Version 1.0 will be released
after the basic analyses (DC/AC/Transient/HB) are implemented and
validated and program interfaces become stable enough.
