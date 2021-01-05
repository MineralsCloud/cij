Output parameters
=================

Use the following keywords to prepare the ``output`` section of your program
configuration file.

.. jinja:: writer_rules

    .. list-table:: Output parameter definitions
        :header-rows: 1
        :stub-columns: 1

        * - Description
          - keywords
          - default unit
          - output filename pattern
    {% for rule in rules %}
        * - {{rule.description}}
          - {% for kw in rule["keywords"] %}``{{kw}}``{{ ", " if not loop.last }}{% endfor %}
          - {{rule.unit}}
          - ``{{rule.fname_pattern}}``
    {% endfor %}

For actually output file name the ``{ij}`` in the output file name pattern will
be replaced by the subscripts of elastic constants symbols in Voigt notation,
and the ``{base}`` will be replaced by ``tv`` and ``tp`` (for value defined on
the grid of :math:`(T,V)` and :math:`(T,P)`, respectively).