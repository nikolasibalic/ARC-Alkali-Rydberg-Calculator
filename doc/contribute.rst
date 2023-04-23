How to contribute
=================

Both new data and calculations for the core library are welcome. We will also
include codes that solve specific research questions in the `arc.advanced <./advanced.html>`_ .
Ideally, this package/library will grow into a community project,
as a community-maintained resource for atomic physics community. Full code is
accessible from `GitHub <https://github.com/nikolasibalic/ARC-Alkali-Rydberg-Calculator>`_, so please fork the project, and submit improvements,
additional modules, new features, or just suggest ideas.
For inspiration, we have a list of possible features for developent below.


Ideas for development
---------------------
This is incomplete list of some of the modules that can be added to the library:

* Photoionization
* Collisional cross-sections
* Dressing potentials
* New methods for precise calculation of DME
* Hyperfine structure calculations
* ... (add your own ideas)

Before you start
----------------

You can contribute code to core of the ARC library, or add specific calculation
example to `arc.advanced`. It is highly recommended, especially if you aim
expanding the core of the library, to open an Issue on GitHub ARC page
and discuss idea with main developers to make sure that good structure is
chosen, that no one else is working at the same idea at the same moment,
and that existing ARC calculations are reused where possible. The quick
discussion should enable smoother development and quicker integration of your
contribution to the main package.

Naming conventions
------------------

For the sake of consistency, readability and cross-linking with the written literature, please follow the following for contributions:

* Names and method/subdivision should reflect **structure of knowledge in atomic physics**, NOT low-level implementation structure.

* Names should be sensible to atomic physicists (even if they are not familiar with coding).

* Use long self-descriptive variable names (so that the code is self-documented and readable in itself) and write short comment on functions of code subsections.

* Use `Google style docstrings for code documentation <https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html>`_ (we are using Sphinx Napoleon extension for generating documentation)

* Add references to original papers in comments and docstrings.

* Ideally, `format your code using black <#code-format>`_

Finally, this is the naming convention. of the original package. For consistency, we suggest following the same naming convention.

* Submodules are lower case, separated by underscore. Example::

    import my_module_name

* Classes are named in CamelCase, for example::

    class MyNewClass:
        ...

* Class methods that return a value start with get in their name, and follow camelCase convention, for example::

    def getSumOfTwoNumbers(a,b):
        return a+b

* Class methods don't return a value are named in camelCase, for example::

    def defineBasis():
        ...

Code format
-----------

Please run ``black``  in the repository folder fixing line length to 80 characters::

    black -l 80 .

And then run linter

    flake8 .

and resolve any erros that appear there before submitting code.