# Posición del TC y del RDSR en DermaDose Surface

El panel de Slicer permite leer el TC seleccionado y seleccionar un RDSR para comparar
su posición nominal. La extracción actualiza el panel automáticamente; la exportación
incluye la evidencia en el JSON del STL. No se añaden matrices a los mensajes de consola.

- TC: conserva IJK→RAS real de Slicer y su conversión explícita IJK→LPS, junto con
  PatientPosition, FrameOfReferenceUID y los campos ImageOrientationPatient y
  ImagePositionPatient disponibles en los DICOM originales, sin leer píxeles.
  Para CT Enhanced cuya geometría se encuentra en grupos funcionales, la geometría
  efectiva del volumen sigue siendo la matriz de Slicer; los campos superiores pueden
  estar ausentes. Un NRRD sin sus DICOM no permite inferir HFS/HFP por su matriz.
- RDSR: lee por evento los conceptos DCM 113743, 113744 y 113745, incluidos modificadores
  anidados. Reconoce códigos SCT y los equivalentes SRT/SNM3 de supino, prono y decúbito
  lateral, más cabeza/pies primero. Conserva los códigos fuente y no interpreta CodeMeaning.
  Valores desconocidos, contradictorios, inclinaciones no caracterizadas o cambios entre
  eventos no se sustituyen por una posición global.
- El STL ya se exporta en RAS/mm: no debe volver a rotarse según PatientPosition del TC.

## Límite deliberado: esto NO establece el registro con el equipo

ImagePositionPatient es una posición en el sistema de referencia del TC. Las posiciones
de mesa del RDSR pertenecen al sistema del equipo intervencionista. Sus orígenes no son
intercambiables, ni compartir postura significa compartir FrameOfReferenceUID.

La nueva lectura no cambia mesh_axes_matrix, mesh_origin_mm ni aplica transformaciones
al cálculo existente. Los JSON indican registration_established=false y transform_applied=false.
Para colocar la malla correctamente todavía se necesita localizar el isocentro en el TC
(o puntos anatómicos de correspondencia) y definir los ejes/convenciones calibrados del
sistema de mesa. Los desplazamientos de mesa podrán aplicarse después respecto a ese
registro, con tratamiento por evento si cambia la postura. No se centra la malla por su
caja envolvente ni se asume que el centro del TC sea el isocentro.

Referencias primarias:
- https://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.7.3.html
- https://dicom.nema.org/medical/dicom/current/output/chtml/part16/sect_CID_19.html
- https://dicom.nema.org/medical/dicom/current/output/chtml/part16/sect_CID_20.html
- https://dicom.nema.org/medical/dicom/current/output/chtml/part16/sect_CID_21.html
- https://dicom.nema.org/medical/dicom/current/output/chtml/part16/sect_TID_10003.html

Las pruebas unitarias cubren códigos, posturas, ausencia de datos, conflictos,
separación por evento y conservación de geometría oblicua sin doble rotación.
La integración de la interfaz requiere comprobarse en Slicer real.
