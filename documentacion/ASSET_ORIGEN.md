# Origen de las mallas humanas de demostración

## Fantoma XCAT elegido (`demo --model mesh`)

El archivo `skindose_verifier/data/xcat_wholebody_body.stl` se extrajo, sin modificar sus bytes, del miembro `WholeBody_Body.stl` de `mesh50XCAT_WholeBody.zip` del repositorio [Mesh-based Human Phantom for Simulation](https://github.com/BenAuer2021/Mesh-based-Human-Phantom-for-Simulation), revisión `1f20b67fd7fa700d27a3dc88ee9ece436658eee0`. El ZIP tiene SHA-256 `08712985e3e3aa30262c25a9f6cdebbe1df95e1b6b3da60167ec95b688ca2061`; el STL extraído tiene SHA-256 `2edd0de03344a1081f2a6567bbcca87b965fbbadd822756e05830cd570cdcb3b`. Se conserva la licencia distribuida en el repositorio en `skindose_verifier/data/XCAT_LICENSE.txt`. El README del autor explica el permiso de distribución gratuita concedido por los creadores del XCAT original y solicita citar sus trabajos y el XCAT original en publicaciones.

Es el contorno corporal **Body** del fantoma mesh50_XCAT de varón estadounidense del percentil 50. Contiene **12 148 triángulos**, es cerrado y conectado, y mide aproximadamente 658 × 280 × 1764 mm. No se importan los órganos internos, de modo que el cálculo es exclusivamente sobre la superficie externa; no calcula dosis de órganos ni hace transporte Monte Carlo. El STL usa coordenadas numéricas de tamaño milimétrico (deducidas de sus cotas y del README; STL no codifica unidades). La demostración declara explícitamente `mesh_unit: mm` y aplica:

```text
x_paciente =  X_malla
y_paciente =  Z_malla - 500 mm
z_paciente = -Y_malla + 40 mm
```

La matriz de ejes y la traslación colocan de forma **aproximada** el tórax en el área irradiada y la cara posterior junto a la mesa sintética; su orientación anterior/posterior y el apoyo son supuestos de demostración que requieren comprobación contra datos anatómicos y de equipo antes de cualquier uso experimental. No representan a un paciente concreto. Los resultados guardan el SHA-256 del STL empleado y las coordenadas ya transformadas de cada triángulo.

## Modelo anterior, disponible como alternativa (`demo --model mesh-pyskindose`)

`skindose_verifier/data/pyskindose_adult_female.stl` procede de `pyskindose/phantom_data/adult_female.stl` del paquete [PySkinDose](https://github.com/rvbCMTS/PySkinDose) 25.1.1. SHA-256 del STL original: `f74681f580aa4b5f281d7f6d49eb9888ffd25c25c9b129f355bf01d0b299cd91`; SHA-256 del wheel: `422396a8717215a460b781d686c1aa98e98b7227307f1faec6073f4e23103804`. Se conserva `PYSKINDOSE_LICENSE.txt`. La [documentación de PySkinDose](https://pyskindose.readthedocs.io/en/latest/pyskindose.html) describe el modelo humano como fantoma MakeHuman. Tiene 26 756 triángulos. Su configuración utiliza unidades de cm y la transformación `x=10X`, `y=10Z+350 mm`, `z=-10Y-100 mm`.

Ambos modelos son geometrías de referencia para pruebas sintéticas; la malla de una persona real requiere segmentación, registro, resolución adecuada y validación física antes de interpretar una PSD como dosis clínica.
