// Appearance / behaviour of the viewer. Edit freely – no other file needs to change.
window.TANK_CONFIG = {
  team: 'Alex West Survey Team',
  title: 'Alex West – Aeration Tank',
  subtitle: 'Preliminary model – verify against the signed drawings before setting out',
  defaults: { scaleToSurvey: true, zExag: 2, cutOn: false, view: '3d', ortho: false },
  // one entry per layer id used by build_data.py. color = base colour, visible = initial state
  layers: {
    raft:     { name: 'Raft (top +1.15)',          color: '#b8bec6', visible: true },
    pcc:      { name: 'PCC (blinding)',            color: '#6f7780', visible: true },
    walls:    { name: 'Walls & end walls',         color: '#8d9aa8', visible: true },
    baffles:  { name: 'Baffle walls',              color: '#5b8fbf', visible: true },
    columns:  { name: 'Columns C3 (baffles)',      color: '#c9a66b', visible: true },
    channels: { name: 'Inlet / outlet channels',   color: '#4ea39a', visible: true },
    supports: { name: 'Channel supports',          color: '#d98c5f', visible: true }
  },
  grids:   { color: '#d8334a', visible: true },    // axes 1..10 and A..F (from the drawing)
  control: { color: '#ffd400', visible: true },    // survey control points P01, P21, P36, P49
  decimals: { en: 3, z: 3, ll: 7 }
};
