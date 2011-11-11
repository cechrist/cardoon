
Coding conventions
==================

General guidelines
++++++++++++++++++

In general all code in the project should follow the following
conventions.  Sometimes a few exceptions may be acceptable if it is
not possible/practical to follow these rules.

1. Use only 80 columns. Break lines if needed.

2. Use lowercase first and (perhaps) underscores for function names:
   copy_to()

3. Classes are capitalized: Element, MyClass

4. Use lowercase first and no underscores for public attributes: vt0,
   isNonlinear

5. Use underscore first for private attributes: _privateVar 

6. module names: first letter lowercase and no underscores: circuit,
   paramset


Documentation
+++++++++++++

Keep doctrings for modules up to date. That documentation will be
included in the module reference. The user manual entry for device
models is automatically generated from the doc string of the
``Device`` class. Similarly, for analyses use the doctring from the
``Analysis`` class in the corresponding module.


