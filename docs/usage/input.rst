Building input
==============

The executation the ``cij`` program requires a configuration file and two
input data files. They are described below.

Program configuration file
""""""""""""""""""""""""""

The configuration file in JSON and YAML format controls the behavior of the program: where to look for the input
data, which parameters to calculate, the pressure and temperature range to work
on.

This files contains three sections, ``qha``, ``elast`` and ``output``. The
``qha`` and ``elast`` will be control the behavior of QHA calculation module and
thermalelasticity calculation module; the ``output`` module controls output
files and formats.

For example in YAML format these sections are written nested:

.. code-block:: yaml

   qha:
     input: # the qha::input
     settings:
       # the qha::settings section
   elast:
     input: # The elast::input
     settings:
       # The the elast::settings section
   output:
       # the output section

In the following table, we summarize the functionalities for these sections in 
the configurational file.

.. jinja:: config_schema

    .. list-table:: **The structure of the configuration file**
        :header-rows: 1
        :stub-columns: 1

        * - key
          - type
          - description
    {% for k, v in flatten_properties(schema).items() %}
        * - ``{{k}}``
          - ``{{v.type}}``
          - {{v.title}}
    {% endfor %}

The following two tables decribes the accepted parameters for ``qha::settings``
and ``elast::settings``.

.. jinja:: config_schema

    {% set schema = schema.definitions.qha_settings %}

    .. list-table:: **{{schema.title}}** (``qha::settings``)
        :header-rows: 1
        :stub-columns: 1

        * - key
          - type
          - description
    {% for k, v in flatten_properties(schema).items() %}
        * - ``{{k}}``
          - ``{{v.type}}``
          - {{v.title}}
    {% endfor %}

.. jinja:: config_schema

    {% set schema = schema.definitions.elast_settings %}

    .. list-table:: **{{schema.title}}** (``elast.settings::settings``)
        :header-rows: 1
        :stub-columns: 1

        * - key
          - type
          - description
    {% for k, v in flatten_properties(schema).items() %}
        * - ``{{k}}``
          - ``{{v.type}}``
          - {{v.title}}
    {% endfor %}


Default values are:

.. include:: ../../cij/data/default/settings.yaml
   :literal:


The QHA input data file
"""""""""""""""""""""""
The first input data is similar to the 
`QHA's input data <https://mineralscloud.github.io/qha/tutorials/run.html#how-to-make-input-data>`_.

The static elastic moduli data file
"""""""""""""""""""""""""""""""""""
