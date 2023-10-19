#############
PyGopher
#############

A Python interface to the PGopher spectral simulation program.

*************
Installation
*************

PyGopher supports Python 3.8+. In a new environment and with the code cloned, you can simply run:

.. code-block:: bash

    pip install .


For development dependencies, run ``pip install -e './[dev]'``

*************
Usage
*************

Simplest way to get started is to do a zero-specification RNG simulation:

.. code-block:: python

    >>> from pygopher import PGopher
    >>> pgo = PGopher.from_rng()
    # returns the line list and partition function as dataframes
    >>> linelist, qpart = pgo.simulate()


More flexibility in the documentation will come later, but for the most part,
we keep the same abstraction as PGopher, and so ``Molecule``, ``Simulation``,
``Species`` objects can be freely configured as if using the GUI by passing
kwargs into ``PGopher.from_rng``:

.. code-block:: python

    >>> pgo = PGopher.from_rng(species_kwargs={...}, mixture_kwargs={...})
