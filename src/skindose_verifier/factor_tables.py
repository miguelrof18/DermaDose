"""Bounded, data-independent multilinear interpolation with source-node audit."""
from dataclasses import dataclass
from itertools import product
import math

import numpy as np

from .domain import InputError


@dataclass(frozen=True)
class FactorPair:
    k_bs: np.ndarray
    k_f: np.ndarray
    audit: dict


class RectilinearFactorTable:
    """A complete Cartesian grid of paired factors; never fills missing nodes."""

    def __init__(self, rows, axes):
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

    def evaluate(self, **coordinates):
        if set(coordinates) != set(self.axes):
            raise InputError(f'Se requieren exactamente los ejes {self.axes}.')
        try:
            values = np.broadcast_arrays(*(np.asarray(coordinates[a], dtype=float) for a in self.axes))
        except (ValueError, TypeError) as exc:
            raise InputError('Coordenadas inválidas en tabla de factores.') from exc
        shape = values[0].shape
        if values[0].size == 0:
            raise InputError('Consulta de factores vacía.')
        brackets = []
        for axis, v in zip(self.axes, values):
            nodes = self.nodes[axis]
            if not np.all(np.isfinite(v)) or np.any(v < nodes[0]) or np.any(v > nodes[-1]):
                raise InputError(f'{axis}: fuera del dominio [{nodes[0]}, {nodes[-1]}]; no se extrapola.')
            upper = np.searchsorted(nodes, v.ravel(), side='left')
            lower = np.maximum(upper-1, 0)
            exact = nodes[upper] == v.ravel()
            lower = np.where(exact, upper, lower)
            span = nodes[upper]-nodes[lower]
            weight = np.divide(v.ravel()-nodes[lower], span, out=np.zeros_like(span), where=span != 0)
            brackets.append((nodes[lower], nodes[upper], weight))
        bs, f = np.zeros(values[0].size), np.zeros(values[0].size)
        used = {}; interpolated = np.zeros(values[0].size, dtype=bool)
        for sides in product((0, 1), repeat=len(self.axes)):
            weights = np.ones(values[0].size)
            coords = []
            for side, (lo, hi, w) in zip(sides, brackets):
                weights *= w if side else 1-w
                coords.append(hi if side else lo)
            for i in np.flatnonzero(weights):
                row = self.rows[tuple(c[i] for c in coords)]
                bs[i] += weights[i]*float(row['k_bs'])
                f[i] += weights[i]*float(row['k_f'])
                used[row['row_id']] = row
                interpolated[i] |= weights[i] != 1
        return FactorPair(bs.reshape(shape), f.reshape(shape), {
            'method':'piecewise_multilinear', 'status':'interpolated' if np.any(interpolated) else 'tabulated',
            'extrapolated':False, 'out_of_domain_policy':'error',
            'query_axes':{a: np.unique(v).tolist() for a,v in zip(self.axes, values)},
            'source_nodes':[used[k] for k in sorted(used)],
            'weight_rule':'product over axes of (x1-x)/(x1-x0) or (x-x0)/(x1-x0); exact axes weight 1',
        })
