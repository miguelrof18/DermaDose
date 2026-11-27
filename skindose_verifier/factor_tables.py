"""Bounded, data-independent multilinear interpolation with source-node audit."""
from dataclasses import dataclass
from itertools import product
import math

import numpy as np

from .domain import InputError

# Queries with more distinct values than this are summarised in the audit
# (count, minimum, maximum). A mesh event can query tens of thousands of local
# field sides; listing all of them made each report grow by megabytes without
# adding traceability, because the source nodes are listed separately.
QUERY_AXIS_LISTING_LIMIT = 32


@dataclass(frozen=True)
class FactorPair:
    k_bs: np.ndarray
    k_f: np.ndarray
    audit: dict


def _query_axis_summary(values: np.ndarray):
    unique = np.unique(values)
    if len(unique) <= QUERY_AXIS_LISTING_LIMIT:
        return unique.tolist()
    return {"count": int(len(unique)), "minimum": float(unique[0]), "maximum": float(unique[-1]),
            "listing": "summarised_more_than_%d_values" % QUERY_AXIS_LISTING_LIMIT}


class RectilinearFactorTable:
    """A complete Cartesian grid of paired factors; never fills missing nodes."""

    def __init__(self, rows, axes, *, interpolate_combined=False):
        self.interpolate_combined = interpolate_combined
        self.axes = tuple(axes)
        if not rows or not self.axes:
            raise InputError('Tabla de factores vacía o sin ejes.')
        self.nodes = {a: np.array(sorted({float(r[a]) for r in rows})) for a in self.axes}
        if any(not np.all(np.isfinite(v)) for v in self.nodes.values()):
            raise InputError('Ejes de factores no finitos.')
        self.rows = {}
        for row in rows:
            key = tuple(float(row[a]) for a in self.axes)
            if key in self.rows:
                raise InputError('Nodos duplicados en tabla de factores.')
            if any(not math.isfinite(float(row[k])) or float(row[k]) <= 0 for k in ('k_bs','k_f')):
                raise InputError('Factores no positivos o no finitos.')
            self.rows[key] = dict(row)
        if len(self.rows) != math.prod(len(v) for v in self.nodes.values()):
            raise InputError('Malla incompleta: no se inventan nodos ausentes.')
        # Dense node arrays allow a vectorised evaluation. They are built once
        # from the same rows and are only an index of the stored values.
        shape = tuple(len(self.nodes[a]) for a in self.axes)
        position = {a: {float(x): i for i, x in enumerate(self.nodes[a])} for a in self.axes}
        self._bs = np.empty(shape)
        self._f = np.empty(shape)
        self._keys = np.empty(shape, dtype=object)
        for key, row in self.rows.items():
            index = tuple(position[a][x] for a, x in zip(self.axes, key))
            self._bs[index] = float(row['k_bs'])
            self._f[index] = float(row['k_f'])
            self._keys[index] = key

    def evaluate(self, **coordinates):
        if set(coordinates) != set(self.axes):
            raise InputError(f'Se requieren exactamente los ejes {self.axes}.')
        try:
            values = np.broadcast_arrays(*(np.asarray(coordinates[a], dtype=float) for a in self.axes))
        except (ValueError, TypeError) as exc:
            raise InputError('Coordenadas inválidas en tabla de factores.') from exc
        shape = values[0].shape
        size = values[0].size
        if size == 0:
            raise InputError('Consulta de factores vacía.')
        brackets = []
        for axis, v in zip(self.axes, values):
            nodes = self.nodes[axis]
            if not np.all(np.isfinite(v)) or np.any(v < nodes[0]) or np.any(v > nodes[-1]):
                raise InputError(f'{axis}: fuera del dominio [{nodes[0]}, {nodes[-1]}]; no se extrapola.')
            flat = v.ravel()
            upper = np.searchsorted(nodes, flat, side='left')
            lower = np.maximum(upper-1, 0)
            exact = nodes[upper] == flat
            lower = np.where(exact, upper, lower)
            span = nodes[upper]-nodes[lower]
            weight = np.divide(flat-nodes[lower], span, out=np.zeros_like(span), where=span != 0)
            brackets.append((lower, upper, weight))
        # Same corner order and the same floating-point operation order as the
        # former per-point loop: results are bitwise identical, only faster.
        bs, f = np.zeros(size), np.zeros(size)
        combined = np.zeros(size)
        interpolated = np.zeros(size, dtype=bool)
        used_corners = []
        for sides in product((0, 1), repeat=len(self.axes)):
            weights = np.ones(size)
            index = []
            for side, (lo, hi, w) in zip(sides, brackets):
                weights *= w if side else 1-w
                index.append(hi if side else lo)
            index = tuple(index)
            active = weights != 0
            if not np.any(active):
                continue
            corner_bs, corner_f = self._bs[index], self._f[index]
            bs = np.where(active, bs + weights*corner_bs, bs)
            f = np.where(active, f + weights*corner_f, f)
            combined = np.where(active, combined + weights*corner_bs*corner_f, combined)
            interpolated |= active & (weights != 1)
            used_corners.append(np.ravel_multi_index(tuple(i[active] for i in index), self._bs.shape))
        flat_used = np.unique(np.concatenate(used_corners))
        used = {}
        for flat_index in flat_used:
            row = self.rows[self._keys.flat[int(flat_index)]]
            used[row['row_id']] = row
        if self.interpolate_combined:
            f = combined/bs
        return FactorPair(bs.reshape(shape), f.reshape(shape), {
            'interpolation_quantity':'combined_C' if self.interpolate_combined else 'separate_factors',
            'method':'piecewise_multilinear', 'status':'interpolated' if np.any(interpolated) else 'tabulated',
            'extrapolated':False, 'out_of_domain_policy':'error',
            'query_axes':{a: _query_axis_summary(v) for a,v in zip(self.axes, values)},
            'source_nodes':[used[k] for k in sorted(used)],
            'weight_rule':'product over axes of (x1-x)/(x1-x0) or (x-x0)/(x1-x0); exact axes weight 1',
        })
