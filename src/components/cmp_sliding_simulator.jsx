import React, { useState, useEffect, useCallback, useRef } from 'react';
import { XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, BarChart, Bar } from 'recharts';

const WaferSimulator = () => {
  const [params, setParams] = useState({
    R_wafer: 150,
    R_pad: 400,
    rpm_pad: 31,
    rpm_wafer: 30,
    R_c: 225,
    time_total: 150,
    num_time_steps: 1500,
    oscillation_enabled: true,
    oscillation_amplitude: 40,
    oscillation_speed: 5.0
  });

  const [results, setResults] = useState(null);
  const [isSimulating, setIsSimulating] = useState(false);

  const [customPoints, setCustomPoints] = useState([]);
  const [newPointRadius, setNewPointRadius] = useState(75);
  const [newPointAngle, setNewPointAngle] = useState(0);
  const [animationTime, setAnimationTime] = useState(0);
  const [isPlaying, setIsPlaying] = useState(false);
  const [animationSpeed, setAnimationSpeed] = useState(1.0);
  const [currentView, setCurrentView] = useState('trajectory');
  const [selectedPadRadii, setSelectedPadRadii] = useState([200]);
  const [barAnalysis, setBarAnalysis] = useState({ 
    padBarData: Array.from({ length: 40 }, (_, i) => ({
      xRange: `${-200 + i * 10}~${-190 + i * 10}`,
      count: 0
    })), 
    waferBarData: Array.from({ length: 30 }, (_, i) => ({
      xRange: `${-150 + i * 10}~${-140 + i * 10}`,
      distance: 0
    }))
  });
  const [simulationKey, setSimulationKey] = useState(0);
  const [isInitialLoad, setIsInitialLoad] = useState(true);

  const [visibleGraphs, setVisibleGraphs] = useState({
    padCoordinate: true,
    heatmap: true,
    velocityVector: true,
    pad200mmTrace: true
  });

  const [activeWafers, setActiveWafers] = useState({
    left: true,
    right: true
  });

  const animationRef = useRef();
  const isPlayingRef = useRef(false);
  const lastFrameTimeRef = useRef(0);

  const animate = useCallback((currentTime) => {
    if (!isPlayingRef.current) return;

    if (lastFrameTimeRef.current === 0) {
      lastFrameTimeRef.current = currentTime;
    }

    const deltaTime = currentTime - lastFrameTimeRef.current;
    lastFrameTimeRef.current = currentTime;

    const realTimeSeconds = deltaTime / 1000;
    const simulationTimeIncrement = realTimeSeconds * animationSpeed * 10;

    setAnimationTime(prev => {
      const next = prev + simulationTimeIncrement;
      const maxTime = results?.timeSeriesData?.length - 1 || 0;
      
      if (next >= maxTime) {
        isPlayingRef.current = false;
        setIsPlaying(false);
        if (animationRef.current) {
          cancelAnimationFrame(animationRef.current);
          animationRef.current = null;
        }
        return maxTime;
      }
      
      if (isPlayingRef.current) {
        animationRef.current = requestAnimationFrame(animate);
      }
      return next;
    });
  }, [animationSpeed, results]);

  const runSimulation = useCallback(() => {
    setIsSimulating(true);

    setTimeout(() => {
      const omega_pad = (params.rpm_pad * 2 * Math.PI) / 60;
      const omega_wafer = (params.rpm_wafer * 2 * Math.PI) / 60;
      
      const time_steps = Array.from({length: params.num_time_steps}, (_, i) => 
        (i * params.time_total) / (params.num_time_steps - 1)
      );
      const dt = time_steps[1] - time_steps[0];

      let wafer_x_oscillation, wafer_vx_oscillation;
      if (params.oscillation_enabled) {
        const freq = params.oscillation_speed / (2 * params.oscillation_amplitude);
        wafer_x_oscillation = time_steps.map(t => 
          params.oscillation_amplitude * Math.sin(2 * Math.PI * freq * t)
        );
        wafer_vx_oscillation = time_steps.map(t => 
          (2 * Math.PI * freq) * params.oscillation_amplitude * Math.cos(2 * Math.PI * freq * t)
        );
      } else {
        wafer_x_oscillation = new Array(params.num_time_steps).fill(0);
        wafer_vx_oscillation = new Array(params.num_time_steps).fill(0);
      }

      const wafer_points = [];
      const target_points = 300;
      const grid_size = Math.ceil(Math.sqrt(target_points * 4 / Math.PI));
      const cell_size = (2 * params.R_wafer) / grid_size;
      
      for (let i = 0; i < grid_size; i++) {
        for (let j = 0; j < grid_size; j++) {
          const x = (i + 0.5) * cell_size - params.R_wafer;
          const y = (j + 0.5) * cell_size - params.R_wafer;
          
          if (Math.sqrt(x * x + y * y) <= params.R_wafer) {
            const r = Math.sqrt(x * x + y * y);
            const theta = Math.atan2(y, x);
            
            let zone = 'outer';
            if (r < params.R_wafer * 0.3) zone = 'inner';
            else if (r < params.R_wafer * 0.7) zone = 'middle';
            
            wafer_points.push({x, y, r, theta, zone});
          }
        }
      }

      // 속도벡터용 점들 - 중심점과 반경 150mm 지점 고정
      const vectorPoints = [];
      
      vectorPoints.push({x: 0, y: 0, r: 0, theta: 0, zone: 'center', vectorId: 0});
      vectorPoints.push({x: 150, y: 0, r: 150, theta: 0, zone: 'edge', vectorId: 1});
      
      const rings = 4;
      const totalVectorPoints = 18;
      
      for (let ring = 1; ring <= rings; ring++) {
        const radius = (ring * params.R_wafer) / rings;
        const pointsInRing = Math.ceil(totalVectorPoints * ring / (rings * (rings + 1) / 2));
        
        for (let i = 0; i < pointsInRing; i++) {
          const angle = (2 * Math.PI * i) / pointsInRing;
          const x = radius * Math.cos(angle);
          const y = radius * Math.sin(angle);
          
          if (Math.sqrt(x * x + y * y) <= params.R_wafer) {
            vectorPoints.push({
              x: x,
              y: y,
              r: radius,
              theta: angle,
              zone: 'vector',
              vectorId: vectorPoints.length
            });
          }
        }
      }
      
      vectorPoints.forEach(point => {
        wafer_points.push(point);
      });

      customPoints.forEach(point => {
        const x = point.radius * Math.cos(point.angle * Math.PI / 180);
        const y = point.radius * Math.sin(point.angle * Math.PI / 180);
        if (Math.sqrt(x * x + y * y) <= params.R_wafer) {
          wafer_points.push({
            x: x, 
            y: y, 
            r: point.radius, 
            theta: point.angle * Math.PI / 180, 
            zone: 'custom',
            customId: point.id
          });
        }
      });

      const num_actual_grid = wafer_points.length;
      const timeSeriesData = [];
      const trajectoryData = wafer_points.map(() => []);
      const cumulativeDistance = wafer_points.map(() => 0);
      const padTrajectoryData = wafer_points.map(() => []);
      const padTrajectoryDataLeft = wafer_points.map(() => []);
      const velocityData = wafer_points.map(() => []);
      
      const padRadiiTrajectories = {};
      selectedPadRadii.forEach(radius => {
        padRadiiTrajectories[radius] = [];
      });
      
      let rotationTrackingByRadius = {};
      selectedPadRadii.forEach(radius => {
        rotationTrackingByRadius[radius] = {
          lastPadAngle: 0,
          currentRotationTraces: [],
          rotationCount: 0,
          allRotations: []
        };
      });

      for (let k = 0; k < params.num_time_steps; k++) {
        const t = time_steps[k];
        const pad_theta = omega_pad * t;
        const wafer_theta = omega_wafer * t;
        
        const wafer_xc = params.R_c + wafer_x_oscillation[k];
        const wafer_yc = 0.0;
        
        const wafer_xc_left = -(params.R_c - wafer_x_oscillation[k]);
        const wafer_yc_left = 0.0;
        
        const cos_w = Math.cos(wafer_theta);
        const sin_w = Math.sin(wafer_theta);

        const timeStepData = {
          time: t,
          wafer_center_x: wafer_xc,
          wafer_center_y: wafer_yc,
          wafer_center_x_left: wafer_xc_left,
          wafer_center_y_left: wafer_yc_left,
          wafer_angle: wafer_theta,
          pad_angle: pad_theta,
          points: []
        };

        selectedPadRadii.forEach(padRadius => {
          const padPoint = { x: padRadius, y: 0 };
          const cos_p = Math.cos(pad_theta);
          const sin_p = Math.sin(pad_theta);
          const global_x = padPoint.x * cos_p - padPoint.y * sin_p;
          const global_y = padPoint.x * sin_p + padPoint.y * cos_p;
          
          const rel_x = global_x - wafer_xc;
          const rel_y = global_y - wafer_yc;
          
          const cos_w_inv = Math.cos(-wafer_theta);
          const sin_w_inv = Math.sin(-wafer_theta);
          const wafer_x = rel_x * cos_w_inv - rel_y * sin_w_inv;
          const wafer_y = rel_x * sin_w_inv + rel_y * cos_w_inv;
          
          const fullRotation = 2 * Math.PI;
          const currentRotation = Math.floor(pad_theta / fullRotation);
          const tracking = rotationTrackingByRadius[padRadius];
          
          if (currentRotation > tracking.rotationCount) {
            if (tracking.currentRotationTraces.length > 0) {
              tracking.allRotations.push({
                rotation: tracking.rotationCount,
                traces: [...tracking.currentRotationTraces]
              });
            }
            tracking.currentRotationTraces = [];
            tracking.rotationCount = currentRotation;
          }
          
          // 웨이퍼 범위 검사를 더 관대하게 조정
          if (Math.sqrt(wafer_x * wafer_x + wafer_y * wafer_y) <= params.R_wafer * 1.1) {
            // 적응형 샘플링: 곡률 변화가 큰 구간에서 더 많은 점 생성
            let shouldAddPoint = tracking.currentRotationTraces.length === 0;
            
            if (!shouldAddPoint && tracking.currentRotationTraces.length > 0) {
              const prevPoint = tracking.currentRotationTraces[tracking.currentRotationTraces.length - 1];
              const distance = Math.sqrt((wafer_x - prevPoint.x)**2 + (wafer_y - prevPoint.y)**2);
              const angle_change = Math.abs(pad_theta - prevPoint.pad_angle);
              
              // 거리 변화 2mm 이상 또는 각도 변화 0.015 라디안 이상시 점 추가
              shouldAddPoint = distance > 2 || angle_change > 0.015;
            }
            
            if (shouldAddPoint) {
              tracking.currentRotationTraces.push({ 
                x: wafer_x, 
                y: wafer_y, 
                time: t, 
                pad_angle: pad_theta,
                step: k,
                radius: padRadius
              });
            }
          }
        });

        for (let i = 0; i < num_actual_grid; i++) {
          const wx = wafer_points[i].x;
          const wy = wafer_points[i].y;
          const zone = wafer_points[i].zone;
          
          const abs_x = wafer_xc + (wx * cos_w - wy * sin_w);
          const abs_y = wafer_yc + (wx * sin_w + wy * cos_w);
          
          const v_center_x = wafer_vx_oscillation[k];
          const v_center_y = 0.0;
          const v_local_x = -omega_wafer * wx * sin_w - omega_wafer * wy * cos_w;
          const v_local_y = omega_wafer * wx * cos_w - omega_wafer * wy * sin_w;
          const v_gx = v_center_x + v_local_x;
          const v_gy = v_center_y + v_local_y;
          
          const v_pad_x = -omega_pad * abs_y;
          const v_pad_y = omega_pad * abs_x;
          
          const v_rel_x = v_gx - v_pad_x;
          const v_rel_y = v_gy - v_pad_y;
          const speed_rel = Math.sqrt(v_rel_x * v_rel_x + v_rel_y * v_rel_y);
          
          if (k > 0) {
            cumulativeDistance[i] += speed_rel * dt;
          }
          
          const cos_p_neg = Math.cos(-pad_theta);
          const sin_p_neg = Math.sin(-pad_theta);
          const pad_local_x = abs_x * cos_p_neg - abs_y * sin_p_neg;
          const pad_local_y = abs_x * sin_p_neg + abs_y * cos_p_neg;
          
          padTrajectoryData[i].push({
            time: t,
            x: pad_local_x,
            y: pad_local_y,
            abs_x: abs_x,
            abs_y: abs_y
          });

          const abs_x_left = wafer_xc_left + (wx * cos_w - wy * sin_w);
          const abs_y_left = wafer_yc_left + (wx * sin_w + wy * cos_w);
          
          const pad_local_x_left = abs_x_left * cos_p_neg - abs_y_left * sin_p_neg;
          const pad_local_y_left = abs_x_left * sin_p_neg + abs_y_left * cos_p_neg;
          
          padTrajectoryDataLeft[i].push({
            time: t,
            x: pad_local_x_left,
            y: pad_local_y_left,
            abs_x: abs_x_left,
            abs_y: abs_y_left
          });
          
          trajectoryData[i].push({
            time: t,
            x: abs_x,
            y: abs_y,
            speed: speed_rel,
            cumulative_distance: cumulativeDistance[i],
            zone: zone
          });

          velocityData[i].push({
            time: t,
            vx: v_rel_x,
            vy: v_rel_y,
            speed: speed_rel
          });
          
          timeStepData.points.push({
            id: i,
            x: abs_x,
            y: abs_y,
            speed: speed_rel,
            cumulative_distance: cumulativeDistance[i],
            zone: zone
          });
        }
        
        timeSeriesData.push(timeStepData);
      }
      
      selectedPadRadii.forEach(padRadius => {
        const tracking = rotationTrackingByRadius[padRadius];
        if (tracking.currentRotationTraces.length > 0) {
          tracking.allRotations.push({
            rotation: tracking.rotationCount,
            traces: [...tracking.currentRotationTraces]
          });
        }
        padRadiiTrajectories[padRadius] = tracking.allRotations;
      });

      const maxDistance = Math.max(...cumulativeDistance);
      const minDistance = Math.min(...cumulativeDistance);
      const avgDistance = cumulativeDistance.reduce((a, b) => a + b, 0) / cumulativeDistance.length;

      setResults({
        wafer_points,
        timeSeriesData,
        trajectoryData,
        padTrajectoryData,
        padTrajectoryDataLeft,
        velocityData,
        padRadiiTrajectories,
        finalDistances: cumulativeDistance,
        statistics: { maxDistance, minDistance, avgDistance },
        R_wafer: params.R_wafer,
        R_pad: params.R_pad,
        time_steps
      });

      setIsSimulating(false);
      setIsInitialLoad(false);
    }, 100);
  }, [params, customPoints, selectedPadRadii]);

  useEffect(() => {
    runSimulation();
  }, [runSimulation, simulationKey]);

  // 📌 실시간 막대그래프 데이터 생성 (애니메이션 시간에 따라 업데이트)
  useEffect(() => {
    if (!results || !results.padTrajectoryData || !results.trajectoryData) return;

    // 현재 애니메이션 시간까지의 데이터만 사용
    const currentTimeIndex = Math.min(Math.floor(animationTime), results.timeSeriesData.length - 1);
    
    // 1. Pad 좌표계 자취 분석 (-200 ~ 200, 10mm 단위)
    const padBins = Array.from({ length: 40 }, (_, i) => -200 + i * 10);
    const padCounts = new Array(padBins.length).fill(0);

    results.padTrajectoryData.forEach((trajectory) => {
      // 현재 시간까지의 궤적만 사용
      const currentTrajectory = trajectory.slice(0, currentTimeIndex + 1);
      currentTrajectory.forEach(({ x }) => {
        const binIdx = Math.floor((x + 200) / 10);
        if (binIdx >= 0 && binIdx < padBins.length) padCounts[binIdx]++;
      });
    });

    const padBarData = padBins.map((x, i) => ({
      xRange: `${x}~${x + 10}`,
      count: padCounts[i]
    }));

    // 2. 웨이퍼 자취 누적 이동거리 분석 (-150 ~ 150)
    const waferBins = Array.from({ length: 30 }, (_, i) => -150 + i * 10);
    const waferDistances = new Array(waferBins.length).fill(0);

    results.trajectoryData.forEach((trajectory, pointIdx) => {
      // 현재 시간까지의 궤적만 사용
      const currentTrajectory = trajectory.slice(0, currentTimeIndex + 1);
      
      // 각 구간별로 이동거리를 계산
      for (let i = 1; i < currentTrajectory.length; i++) {
        const prevPoint = currentTrajectory[i - 1];
        const currPoint = currentTrajectory[i];
        
        // 이전 점과 현재 점의 X 좌표
        const prevX = prevPoint.x;
        const currX = currPoint.x;
        
        // 이동거리 계산
        const distance = Math.sqrt(
          Math.pow(currPoint.x - prevPoint.x, 2) + 
          Math.pow(currPoint.y - prevPoint.y, 2)
        );
        
        // 이전 점과 현재 점이 같은 구간에 있는지 확인
        const prevBinIdx = Math.floor((prevX + 150) / 10);
        const currBinIdx = Math.floor((currX + 150) / 10);
        
        if (prevBinIdx === currBinIdx && prevBinIdx >= 0 && prevBinIdx < waferBins.length) {
          // 같은 구간 내에서의 이동
          waferDistances[prevBinIdx] += distance;
        } else {
          // 구간을 건너는 이동의 경우, 각 구간에 비례하여 분배
          if (prevBinIdx >= 0 && prevBinIdx < waferBins.length) {
            waferDistances[prevBinIdx] += distance * 0.5;
          }
          if (currBinIdx >= 0 && currBinIdx < waferBins.length) {
            waferDistances[currBinIdx] += distance * 0.5;
          }
        }
      }
    });

    const waferBarData = waferBins.map((x, i) => ({
      xRange: `${x}~${x + 10}`,
      distance: waferDistances[i]
    }));

    setBarAnalysis({ padBarData, waferBarData });
  }, [results, animationTime, isInitialLoad]);

  const addCustomPoint = () => {
    if (newPointRadius <= params.R_wafer && newPointRadius >= 0) {
      const newPoint = {
        id: Date.now(),
        radius: newPointRadius,
        angle: newPointAngle
      };
      setCustomPoints(prev => [...prev, newPoint]);
    }
  };

  const removeCustomPoint = (pointId) => {
    setCustomPoints(prev => prev.filter(p => p.id !== pointId));
  };

  const toggleAnimation = useCallback(() => {
    if (isPlaying) {
      isPlayingRef.current = false;
      setIsPlaying(false);
      if (animationRef.current) {
        cancelAnimationFrame(animationRef.current);
        animationRef.current = null;
      }
    } else {
      if (!results || !results.timeSeriesData || results.timeSeriesData.length === 0) return;
      isPlayingRef.current = true;
      setIsPlaying(true);
    }
  }, [isPlaying, results]);

  useEffect(() => {
    if (!isPlaying || !results || !results.timeSeriesData || results.timeSeriesData.length === 0) {
      isPlayingRef.current = false;
      if (animationRef.current) {
        cancelAnimationFrame(animationRef.current);
        animationRef.current = null;
      }
      return;
    }

    if (isPlaying && isPlayingRef.current) {
      lastFrameTimeRef.current = 0;
      animationRef.current = requestAnimationFrame(animate);
    }

    return () => {
      isPlayingRef.current = false;
      if (animationRef.current) {
        cancelAnimationFrame(animationRef.current);
        animationRef.current = null;
      }
    };
  }, [isPlaying, animate]);

  const resetAnimation = useCallback(() => {
    isPlayingRef.current = false;
    setIsPlaying(false);
    if (animationRef.current) {
      cancelAnimationFrame(animationRef.current);
      animationRef.current = null;
    }
    lastFrameTimeRef.current = 0;
    setAnimationTime(0);
    // 막대그래프 데이터 초기화
    setBarAnalysis({ 
      padBarData: Array.from({ length: 40 }, (_, i) => ({
        xRange: `${-200 + i * 10}~${-190 + i * 10}`,
        count: 0
      })), 
      waferBarData: Array.from({ length: 30 }, (_, i) => ({
        xRange: `${-150 + i * 10}~${-140 + i * 10}`,
        distance: 0
      }))
    });
    // 초기 로딩 상태로 재설정
    setIsInitialLoad(true);
    // 시뮬레이션 키 증가하여 강제 재실행
    setSimulationKey(prev => prev + 1);
  }, []);

  const toggleGraph = (graphType) => {
    setVisibleGraphs(prev => ({
      ...prev,
      [graphType]: !prev[graphType]
    }));
  };

  if (!results || !results.timeSeriesData || results.timeSeriesData.length === 0) {
    return (
      <div className="max-w-7xl mx-auto p-4 sm:p-6 space-y-6">
        <h1 className="text-2xl sm:text-3xl font-bold text-center">웨이퍼 폴리싱 시뮬레이터</h1>
        
        <div className="bg-white p-4 rounded-lg shadow-lg space-y-4">
          <h2 className="text-lg font-bold">시뮬레이션 제어</h2>
          <div className="grid grid-cols-2 sm:grid-cols-4 gap-2">
            <div>
              <label className="block text-xs font-medium mb-1">웨이퍼R</label>
              <input
                type="number"
                value={params.R_wafer}
                onChange={(e) => setParams({...params, R_wafer: Number(e.target.value)})}
                className="w-full px-2 py-1 text-sm border rounded"
              />
            </div>
            <div>
              <label className="block text-xs font-medium mb-1">패드R</label>
              <input
                type="number"
                value={params.R_pad}
                onChange={(e) => setParams({...params, R_pad: Number(e.target.value)})}
                className="w-full px-2 py-1 text-sm border rounded"
              />
            </div>
            <div>
              <label className="block text-xs font-medium mb-1">패드RPM</label>
              <input
                type="number"
                value={params.rpm_pad}
                onChange={(e) => setParams({...params, rpm_pad: Number(e.target.value)})}
                className="w-full px-2 py-1 text-sm border rounded"
              />
            </div>
            <div>
              <label className="block text-xs font-medium mb-1">웨이퍼RPM</label>
              <input
                type="number"
                value={params.rpm_wafer}
                onChange={(e) => setParams({...params, rpm_wafer: Number(e.target.value)})}
                className="w-full px-2 py-1 text-sm border rounded"
              />
            </div>
          </div>
        </div>

        {isSimulating && (
          <div className="w-full h-96 bg-gray-100 border rounded-lg flex items-center justify-center">
            <div className="text-center">
              <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-blue-500 mx-auto mb-4"></div>
              <p className="text-gray-600">시뮬레이션 실행 중...</p>
            </div>
          </div>
        )}
      </div>
    );
  }

  const currentTimeIndex = Math.min(Math.floor(animationTime), results.timeSeriesData.length - 1);
  const currentData = results.timeSeriesData[currentTimeIndex];

  return (
    <div className="max-w-7xl mx-auto p-4 sm:p-6 space-y-6">
      <h1 className="text-2xl sm:text-3xl font-bold text-center">웨이퍼 폴리싱 시뮬레이터</h1>

      <div className="bg-white p-4 rounded-lg shadow-lg space-y-4">
        <h2 className="text-lg font-bold">시뮬레이션 제어</h2>

        <div className="grid grid-cols-2 sm:grid-cols-4 gap-2">
          <div>
            <label className="block text-xs font-medium mb-1">웨이퍼R</label>
            <input
              type="number"
              value={params.R_wafer}
              onChange={(e) => setParams({...params, R_wafer: Number(e.target.value)})}
              className="w-full px-2 py-1 text-sm border rounded"
            />
          </div>
          <div>
            <label className="block text-xs font-medium mb-1">패드R</label>
            <input
              type="number"
              value={params.R_pad}
              onChange={(e) => setParams({...params, R_pad: Number(e.target.value)})}
              className="w-full px-2 py-1 text-sm border rounded"
            />
          </div>
          <div>
            <label className="block text-xs font-medium mb-1">패드RPM</label>
            <input
              type="number"
              value={params.rpm_pad}
              onChange={(e) => setParams({...params, rpm_pad: Number(e.target.value)})}
              className="w-full px-2 py-1 text-sm border rounded"
            />
          </div>
          <div>
            <label className="block text-xs font-medium mb-1">웨이퍼RPM</label>
            <input
              type="number"
              value={params.rpm_wafer}
              onChange={(e) => setParams({...params, rpm_wafer: Number(e.target.value)})}
              className="w-full px-2 py-1 text-sm border rounded"
            />
          </div>
        </div>

        <div className="border-t pt-3">
          <h3 className="font-semibold text-sm mb-2">진동 제어</h3>
          <div className="grid grid-cols-2 sm:grid-cols-3 gap-2 mb-2">
            <div>
              <label className="block text-xs font-medium mb-1">진동폭(mm)</label>
              <input
                type="number"
                value={params.oscillation_amplitude}
                onChange={(e) => setParams({...params, oscillation_amplitude: Number(e.target.value)})}
                className="w-full px-2 py-1 text-sm border rounded"
              />
            </div>
            <div>
              <label className="block text-xs font-medium mb-1">진동속도</label>
              <input
                type="number"
                step="0.1"
                value={params.oscillation_speed}
                onChange={(e) => setParams({...params, oscillation_speed: Number(e.target.value)})}
                className="w-full px-2 py-1 text-sm border rounded"
              />
            </div>
            <div className="flex items-center">
              <label className="flex items-center">
                <input
                  type="checkbox"
                  checked={params.oscillation_enabled}
                  onChange={(e) => setParams({...params, oscillation_enabled: e.target.checked})}
                  className="mr-2"
                />
                진동 활성화
              </label>
            </div>
          </div>
        </div>

        <div className="border-t pt-3">
          <h3 className="font-semibold text-sm mb-2">웨이퍼 선택</h3>
          <div className="flex items-center space-x-4 mb-3">
            <label className="flex items-center">
              <input
                type="checkbox"
                checked={activeWafers.left}
                onChange={() => setActiveWafers(prev => ({...prev, left: !prev.left}))}
                className="mr-2"
              />
              L 웨이퍼 (좌측/녹색)
            </label>
            <label className="flex items-center">
              <input
                type="checkbox"
                checked={activeWafers.right}
                onChange={() => setActiveWafers(prev => ({...prev, right: !prev.right}))}
                className="mr-2"
              />
              R 웨이퍼 (우측/파란색)
            </label>
          </div>
        </div>

        <div className="border-t pt-3">
          <h3 className="font-semibold text-sm mb-2">그래프 선택</h3>
          <div className="grid grid-cols-2 gap-2">
            <label className="flex items-center">
              <input
                type="checkbox"
                checked={visibleGraphs.padCoordinate}
                onChange={() => toggleGraph('padCoordinate')}
                className="mr-2"
              />
              패드 좌표계
            </label>
            <label className="flex items-center">
              <input
                type="checkbox"
                checked={visibleGraphs.heatmap}
                onChange={() => toggleGraph('heatmap')}
                className="mr-2"
              />
              히트맵
            </label>
            <label className="flex items-center">
              <input
                type="checkbox"
                checked={visibleGraphs.velocityVector}
                onChange={() => toggleGraph('velocityVector')}
                className="mr-2"
              />
              속도벡터
            </label>
            <label className="flex items-center">
              <input
                type="checkbox"
                checked={visibleGraphs.pad200mmTrace}
                onChange={() => toggleGraph('pad200mmTrace')}
                className="mr-2"
              />
              패드자취
            </label>
          </div>
        </div>

        <div className="border-t pt-3">
          <h3 className="font-semibold text-sm mb-2">패드 반지름 선택 (자취용)</h3>
          <div className="grid grid-cols-4 gap-2 mb-2">
            {[50, 100, 150, 200, 250, 300, 350, 390].map(radius => (
              <label key={radius} className="flex items-center">
                <input
                  type="checkbox"
                  checked={selectedPadRadii.includes(radius)}
                  onChange={() => {
                    setSelectedPadRadii(prev => 
                      prev.includes(radius) 
                        ? prev.filter(r => r !== radius)
                        : [...prev, radius].sort((a, b) => a - b)
                    );
                  }}
                  className="mr-1"
                />
                <span className="text-xs">{radius}mm</span>
              </label>
            ))}
          </div>
        </div>

        <div className="border-t pt-3">
          <h3 className="font-semibold text-sm mb-2">추적점 추가</h3>
          <div className="flex items-center space-x-2 mb-2">
            <input
              type="number"
              placeholder="반지름"
              value={newPointRadius}
              onChange={(e) => setNewPointRadius(Number(e.target.value))}
              className="w-20 px-2 py-1 text-sm border rounded"
            />
            <input
              type="number"
              placeholder="각도"
              value={newPointAngle}
              onChange={(e) => setNewPointAngle(Number(e.target.value))}
              className="w-20 px-2 py-1 text-sm border rounded"
            />
            <button
              onClick={addCustomPoint}
              className="px-3 py-1 bg-blue-500 text-white text-sm rounded"
            >
              추가
            </button>
          </div>
          
          {customPoints.length > 0 && (
            <div className="space-y-1">
              {customPoints.map(point => (
                <div key={point.id} className="flex items-center justify-between bg-gray-50 px-2 py-1 rounded text-xs">
                  <span>R: {point.radius}mm, θ: {point.angle}°</span>
                  <button
                    onClick={() => removeCustomPoint(point.id)}
                    className="text-red-500 hover:text-red-700"
                  >
                    삭제
                  </button>
                </div>
              ))}
            </div>
          )}
        </div>
      </div>
      
      <div>
        <div className="flex gap-2 mb-4">
          <button
            onClick={() => setCurrentView('trajectory')}
            className={`px-3 py-2 rounded text-sm ${currentView === 'trajectory' ? 'bg-blue-500 text-white' : 'bg-gray-200'}`}
          >
            실시간 애니메이션
          </button>

        </div>

        {currentView === 'trajectory' && (
          <div className="space-y-6">
            <div className="flex flex-col sm:flex-row justify-between items-start sm:items-center space-y-4 sm:space-y-0">
              <h3 className="text-xl font-semibold">실시간 시뮬레이션</h3>
              <div className="flex flex-col sm:flex-row items-start sm:items-center space-y-2 sm:space-y-0 sm:space-x-4">
                <div className="flex items-center space-x-2">
                  <button
                    onClick={toggleAnimation}
                    className={`px-4 py-2 rounded ${isPlaying ? 'bg-red-500 text-white' : 'bg-green-500 text-white'}`}
                  >
                    {isPlaying ? '일시정지' : '재생'}
                  </button>
                  <button
                    onClick={resetAnimation}
                    className="px-4 py-2 rounded bg-gray-500 text-white"
                  >
                    리셋
                  </button>
                </div>
                <div className="flex items-center space-x-2">
                  <label className="text-sm">속도:</label>
                  <select
                    value={animationSpeed}
                    onChange={(e) => setAnimationSpeed(Number(e.target.value))}
                    className="px-2 py-1 border rounded text-sm"
                  >
                    <option value={0.25}>0.25x</option>
                    <option value={0.5}>0.5x</option>
                    <option value={1.0}>1.0x</option>
                    <option value={2.0}>2.0x</option>
                    <option value={5.0}>5.0x</option>
                  </select>
                </div>
                <span className="text-sm">시간: {currentData.time.toFixed(2)}s</span>
              </div>
            </div>

            <div className="space-y-6">
              {visibleGraphs.padCoordinate && (
                <div className="grid grid-cols-1 lg:grid-cols-2 gap-4">
                  <div className="bg-white border rounded-lg p-4">
                    <h4 className="text-lg font-semibold mb-4">패드 좌표계</h4>
                    <div className="relative w-full h-80 bg-gray-50 border">
                      <svg width="100%" height="100%" viewBox={`${-results.R_pad*1.2} ${-results.R_pad*1.2} ${results.R_pad*2.4} ${results.R_pad*2.4}`}>
                        <circle cx="0" cy="0" r={results.R_pad} fill="lightgray" opacity="0.3" stroke="gray" strokeWidth="2" />
                        <circle cx="0" cy="0" r={params.R_c} fill="none" stroke="orange" strokeWidth="2" strokeDasharray="5,5" opacity="0.7" />
                        
                        <g transform={`rotate(${currentData.pad_angle * 180 / Math.PI})`}>
                          <line x1="0" y1="0" x2={results.R_pad * 0.9} y2="0" stroke="red" strokeWidth="4" opacity="0.7" />
                          <circle cx={results.R_pad * 0.85} cy="0" r="8" fill="red" opacity="0.8" />
                        </g>
                        
                        {activeWafers.right && (
                          <g transform={`translate(${currentData.wafer_center_x}, ${currentData.wafer_center_y})`}>
                            <circle cx="0" cy="0" r={results.R_wafer} fill="none" stroke="blue" strokeWidth="2" opacity="0.4" />
                            <g transform={`rotate(${currentData.wafer_angle * 180 / Math.PI})`}>
                              <line x1="0" y1="0" x2={results.R_wafer * 0.8} y2="0" stroke="blue" strokeWidth="3" />
                            </g>
                            <text x="0" y={results.R_wafer + 20} textAnchor="middle" fontSize="10" fill="blue">R웨이퍼</text>
                          </g>
                        )}

                        {activeWafers.left && (
                          <g transform={`translate(${currentData.wafer_center_x_left}, ${currentData.wafer_center_y_left})`}>
                            <circle cx="0" cy="0" r={results.R_wafer} fill="none" stroke="green" strokeWidth="2" opacity="0.4" />
                            <g transform={`rotate(${currentData.wafer_angle * 180 / Math.PI})`}>
                              <line x1="0" y1="0" x2={results.R_wafer * 0.8} y2="0" stroke="green" strokeWidth="3" />
                            </g>
                            <text x="0" y={results.R_wafer + 20} textAnchor="middle" fontSize="10" fill="green">L웨이퍼</text>
                          </g>
                        )}

                        {activeWafers.right && customPoints.map(point => {
                          const pointIndex = results.wafer_points.findIndex(p => p.customId === point.id);
                          if (pointIndex === -1) return null;
                          const trajectory = results.padTrajectoryData[pointIndex]
                            .slice(0, currentTimeIndex + 1)
                            .filter(trajPoint => Math.sqrt(trajPoint.x * trajPoint.x + trajPoint.y * trajPoint.y) <= results.R_pad);
                          if (trajectory.length < 2) return null;
                          const pathData = trajectory.map((seg, i) => `${i === 0 ? 'M' : 'L'} ${seg.x} ${-seg.y}`).join(' ');
                          return (
                            <g key={`trail-right-${point.id}`}>
                              <path d={pathData} fill="none" stroke="blue" strokeWidth="2" opacity="0.8" />
                              {trajectory.length > 0 && (
                                <circle cx={trajectory[trajectory.length - 1].x} cy={-trajectory[trajectory.length - 1].y} r="4" fill="blue" />
                              )}
                            </g>
                          );
                        })}

                        {activeWafers.left && results.padTrajectoryDataLeft && customPoints.map(point => {
                          const pointIndex = results.wafer_points.findIndex(p => p.customId === point.id);
                          if (pointIndex === -1 || !results.padTrajectoryDataLeft[pointIndex]) return null;
                          const trajectory = results.padTrajectoryDataLeft[pointIndex]
                            .slice(0, currentTimeIndex + 1)
                            .filter(trajPoint => Math.sqrt(trajPoint.x * trajPoint.x + trajPoint.y * trajPoint.y) <= results.R_pad);
                          if (trajectory.length < 2) return null;
                          const pathData = trajectory.map((seg, i) => `${i === 0 ? 'M' : 'L'} ${seg.x} ${-seg.y}`).join(' ');
                          return (
                            <g key={`trail-left-${point.id}`}>
                              <path d={pathData} fill="none" stroke="green" strokeWidth="2" opacity="0.8" />
                              {trajectory.length > 0 && (
                                <circle cx={trajectory[trajectory.length - 1].x} cy={-trajectory[trajectory.length - 1].y} r="4" fill="green" />
                              )}
                            </g>
                          );
                        })}
                      </svg>
                    </div>
                  </div>

                  {/* ✅ 패드 좌표계 X축 누적횟수 그래프 */}
                  <div className="bg-white border rounded-lg p-4">
                    <h4 className="text-lg font-semibold mb-4">Pad profile</h4>
                    <ResponsiveContainer width="100%" height={320}>
                      <BarChart data={barAnalysis.padBarData}>
                        <CartesianGrid strokeDasharray="3 3" />
                        <XAxis dataKey="xRange" angle={-45} textAnchor="end" height={60} />
                        <YAxis />
                        <Tooltip />
                        <Bar dataKey="count" fill="#4ECDC4" />
                      </BarChart>
                    </ResponsiveContainer>
                  </div>
                </div>
              )}
              
              {visibleGraphs.heatmap && (
                <div className="grid grid-cols-1 lg:grid-cols-2 gap-4">
                  <div className="bg-white border rounded-lg p-4">
                    <h4 className="text-lg font-semibold mb-4">Wafer Map</h4>
                    <div className="relative w-full h-80 bg-gray-50 border">
                      <svg width="100%" height="100%" viewBox={`${-results.R_wafer*1.2} ${-results.R_wafer*1.2} ${results.R_wafer*2.4} ${results.R_wafer*2.4}`}>
                        <circle cx="0" cy="0" r={results.R_wafer} fill="none" stroke="blue" strokeWidth="3" />
                        {results.wafer_points.map((point, i) => {
                          const currentDistance = currentTimeIndex < results.trajectoryData[i].length ? 
                            results.trajectoryData[i][currentTimeIndex].cumulative_distance : 0;
                          const allCurrentDistances = results.wafer_points.map((_, idx) => 
                            currentTimeIndex < results.trajectoryData[idx].length ? 
                            results.trajectoryData[idx][currentTimeIndex].cumulative_distance : 0
                          );
                          const maxCurrentDistance = Math.max(...allCurrentDistances);
                          const minCurrentDistance = Math.min(...allCurrentDistances);
                          const range = maxCurrentDistance - minCurrentDistance;
                          const intensity = range > 0 ? (currentDistance - minCurrentDistance) / range : 0;
                          let red, green, blue;
                          if (intensity < 0.25) {
                            red = 0; green = Math.floor(intensity * 4 * 255); blue = 255;
                          } else if (intensity < 0.5) {
                            red = 0; green = 255; blue = Math.floor(255 - (intensity - 0.25) * 4 * 255);
                          } else if (intensity < 0.75) {
                            red = Math.floor((intensity - 0.5) * 4 * 255); green = 255; blue = 0;
                          } else {
                            red = 255; green = Math.floor(255 - (intensity - 0.75) * 4 * 255); blue = 0;
                          }
                          return (
                            <circle key={i} cx={point.x} cy={-point.y} r="9" fill={`rgb(${red}, ${green}, ${blue})`} opacity="0.85" stroke="rgba(0,0,0,0.2)" strokeWidth="0.5" />
                          );
                        })}
                      </svg>
                    </div>
                  </div>

                  {/* ✅ 웨이퍼 히트맵 X축 누적 이동거리 그래프 */}
                  <div className="bg-white border rounded-lg p-4">
                    <h4 className="text-lg font-semibold mb-4">Wafer Profile</h4>
                    <ResponsiveContainer width="100%" height={320}>
                      <BarChart data={barAnalysis.waferBarData}>
                        <CartesianGrid strokeDasharray="3 3" />
                        <XAxis dataKey="xRange" angle={-45} textAnchor="end" height={60} />
                        <YAxis />
                        <Tooltip />
                        <Bar dataKey="distance" fill="#FF6B6B" />
                      </BarChart>
                    </ResponsiveContainer>
                  </div>
                </div>
              )}
              
              {visibleGraphs.velocityVector && (
                <div className="grid grid-cols-1 lg:grid-cols-2 gap-4">
                  <div className="bg-white border rounded-lg p-4">
                    <h4 className="text-lg font-semibold mb-4">속도 벡터</h4>
                    <div className="relative w-full h-80 bg-gray-50 border">
                      <svg width="100%" height="100%" viewBox={`${-results.R_wafer*1.2} ${-results.R_wafer*1.2} ${results.R_wafer*2.4} ${results.R_wafer*2.4}`}>
                        <circle cx="0" cy="0" r={results.R_wafer} fill="none" stroke="blue" strokeWidth="3" />
                        <defs>
                          <marker id="arrowhead" markerWidth="10" markerHeight="7" refX="9" refY="3.5" orient="auto">
                            <polygon points="0 0, 10 3.5, 0 7" fill="red" />
                          </marker>
                        </defs>

                        {/* 고정 벡터 점들 */}
                        {results.wafer_points.filter(point => point.zone === 'center' || point.zone === 'edge' || point.zone === 'vector').map((point, i) => {
                          const velocityIndex = results.wafer_points.findIndex(p => p.x === point.x && p.y === point.y);
                          if (velocityIndex === -1 || currentTimeIndex >= results.velocityData[velocityIndex].length) return null;
                          
                          const velocity = results.velocityData[velocityIndex][currentTimeIndex];
                          const speed = velocity.speed;
                          const scale = Math.min(speed / 100, 50); // 속도에 따른 벡터 크기 조정
                          
                          if (speed < 0.1) return null; // 너무 작은 속도는 표시하지 않음
                          
                          const endX = point.x + (velocity.vx / speed) * scale;
                          const endY = -point.y - (velocity.vy / speed) * scale;
                          
                          return (
                            <g key={`vector-${i}`}>
                              <line 
                                x1={point.x} 
                                y1={-point.y} 
                                x2={endX} 
                                y2={endY} 
                                stroke="red" 
                                strokeWidth="2" 
                                markerEnd="url(#arrowhead)"
                              />
                              <circle 
                                cx={point.x} 
                                cy={-point.y} 
                                r="3" 
                                fill={point.zone === 'center' ? 'red' : point.zone === 'edge' ? 'blue' : 'green'} 
                                stroke="white" 
                                strokeWidth="1"
                              />
                            </g>
                          );
                        })}

                        {/* 커스텀 점들의 속도 벡터 */}
                        {customPoints.map(point => {
                          const pointIndex = results.wafer_points.findIndex(p => p.customId === point.id);
                          if (pointIndex === -1 || currentTimeIndex >= results.velocityData[pointIndex].length) return null;
                          
                          const velocity = results.velocityData[pointIndex][currentTimeIndex];
                          const speed = velocity.speed;
                          const scale = Math.min(speed / 100, 50);
                          
                          if (speed < 0.1) return null;
                          
                          const endX = point.radius * Math.cos(point.angle * Math.PI / 180) + (velocity.vx / speed) * scale;
                          const endY = -(point.radius * Math.sin(point.angle * Math.PI / 180)) - (velocity.vy / speed) * scale;
                          
                          return (
                            <g key={`custom-vector-${point.id}`}>
                              <line 
                                x1={point.radius * Math.cos(point.angle * Math.PI / 180)} 
                                y1={-(point.radius * Math.sin(point.angle * Math.PI / 180))} 
                                x2={endX} 
                                y2={endY} 
                                stroke="purple" 
                                strokeWidth="2.5" 
                                markerEnd="url(#arrowhead)"
                              />
                              <circle 
                                cx={point.radius * Math.cos(point.angle * Math.PI / 180)} 
                                cy={-(point.radius * Math.sin(point.angle * Math.PI / 180))} 
                                r="4" 
                                fill="purple" 
                                stroke="white" 
                                strokeWidth="1"
                              />
                            </g>
                          );
                        })}

                        {/* 범례 */}
                        <g transform="translate(-400, -350)">
                          <rect x="0" y="0" width="120" height="80" fill="white" opacity="0.95" stroke="gray" strokeWidth="1" rx="3"/>
                          <text x="5" y="15" fontSize="10" fontWeight="bold" fill="black">속도 벡터 범례</text>
                          <g>
                            <line x1="5" y1="25" x2="15" y2="25" stroke="red" strokeWidth="2" markerEnd="url(#arrowhead)"/>
                            <circle cx="5" cy="25" r="2" fill="red"/>
                            <text x="20" y="30" fontSize="8" fill="black">중심점</text>
                          </g>
                          <g>
                            <line x1="5" y1="40" x2="15" y2="40" stroke="red" strokeWidth="2" markerEnd="url(#arrowhead)"/>
                            <circle cx="5" cy="40" r="2" fill="blue"/>
                            <text x="20" y="45" fontSize="8" fill="black">가장자리</text>
                          </g>
                          <g>
                            <line x1="5" y1="55" x2="15" y2="55" stroke="red" strokeWidth="2" markerEnd="url(#arrowhead)"/>
                            <circle cx="5" cy="55" r="2" fill="green"/>
                            <text x="20" y="60" fontSize="8" fill="black">고정점</text>
                          </g>
                          <g>
                            <line x1="5" y1="70" x2="15" y2="70" stroke="purple" strokeWidth="2" markerEnd="url(#arrowhead)"/>
                            <circle cx="5" cy="70" r="2" fill="purple"/>
                            <text x="20" y="75" fontSize="8" fill="black">커스텀점</text>
                          </g>
                        </g>
                      </svg>
                    </div>
                  </div>

                  {/* 속도벡터용 점추가 컨트롤 */}
                  <div className="bg-white border rounded-lg p-4">
                    <h4 className="text-lg font-semibold mb-4">속도벡터 점 관리</h4>
                    <div className="space-y-4">
                      <div className="grid grid-cols-2 gap-2">
                        <div>
                          <label className="block text-xs font-medium mb-1">반지름 (mm)</label>
                          <input
                            type="number"
                            placeholder="0-150"
                            value={newPointRadius}
                            onChange={(e) => setNewPointRadius(Number(e.target.value))}
                            className="w-full px-2 py-1 text-sm border rounded"
                          />
                        </div>
                        <div>
                          <label className="block text-xs font-medium mb-1">각도 (도)</label>
                          <input
                            type="number"
                            placeholder="0-360"
                            value={newPointAngle}
                            onChange={(e) => setNewPointAngle(Number(e.target.value))}
                            className="w-full px-2 py-1 text-sm border rounded"
                          />
                        </div>
                      </div>
                      
                      <button
                        onClick={addCustomPoint}
                        className="w-full px-3 py-2 bg-blue-500 text-white text-sm rounded hover:bg-blue-600"
                      >
                        속도벡터 점 추가
                      </button>
                      
                      {customPoints.length > 0 && (
                        <div className="space-y-2">
                          <h5 className="text-sm font-medium">추가된 점들:</h5>
                          <div className="max-h-32 overflow-y-auto space-y-1">
                            {customPoints.map(point => (
                              <div key={point.id} className="flex items-center justify-between bg-gray-50 px-2 py-1 rounded text-xs">
                                <span>R: {point.radius}mm, θ: {point.angle}°</span>
                                <button
                                  onClick={() => removeCustomPoint(point.id)}
                                  className="text-red-500 hover:text-red-700"
                                >
                                  삭제
                                </button>
                              </div>
                            ))}
                          </div>
                        </div>
                      )}
                    </div>
                  </div>
                </div>
              )}
              
              {visibleGraphs.pad200mmTrace && (
                <div className="bg-white border rounded-lg p-4">
                  <h4 className="text-lg font-semibold mb-4">Scratch Trace</h4>
                  <div className="relative w-full h-80 bg-gray-50 border">
                    <svg width="100%" height="100%" viewBox={`${-results.R_wafer*1.2} ${-results.R_wafer*1.2} ${results.R_wafer*2.4} ${results.R_wafer*2.4}`}>
                      <circle cx="0" cy="0" r={results.R_wafer} fill="none" stroke="blue" strokeWidth="3" />
                      
                      {selectedPadRadii.map((padRadius) => {
                        const trajectoryData = results.padRadiiTrajectories[padRadius] || [];
                        const currentPadAngle = results.timeSeriesData[currentTimeIndex].pad_angle;
                        const currentRotation = Math.floor(currentPadAngle / (2 * Math.PI));
                        
                        const radiusColors = {
                          50: '#FF6B6B', 100: '#4ECDC4', 150: '#45B7D1', 200: '#FFA726',
                          250: '#66BB6A', 300: '#AB47BC', 350: '#FF7043', 390: '#8D6E63'
                        };
                        
                        return trajectoryData.map((rotationData, rotIndex) => {
                          if (rotationData.rotation > currentRotation) return null;
                          
                          let traces = rotationData.traces;
                          if (rotationData.rotation === currentRotation) {
                            traces = traces.filter(point => point.step <= currentTimeIndex);
                          }
                          
                          if (traces.length < 2) return null;
                          
                          // 3차 베지어 곡선으로 부드러운 궤적 생성
                          let pathData = `M ${traces[0].x} ${-traces[0].y}`;
                          
                          if (traces.length === 2) {
                            // 점이 2개일 때는 직선
                            pathData += ` L ${traces[1].x} ${-traces[1].y}`;
                          } else {
                            // 3개 이상일 때는 3차 베지어 곡선
                            for (let i = 1; i < traces.length; i++) {
                              const p0 = i > 0 ? traces[i-1] : traces[i];
                              const p1 = traces[i];
                              const p2 = i < traces.length - 1 ? traces[i+1] : traces[i];
                              const p3 = i < traces.length - 2 ? traces[i+2] : p2;
                              
                              if (i === 1) {
                                // 첫 번째 곡선: 부드러운 시작
                                const cp1x = p0.x + (p1.x - p0.x) * 0.3;
                                const cp1y = -p0.y - (p1.y - p0.y) * 0.3;
                                const cp2x = p1.x - (p2.x - p0.x) * 0.1;
                                const cp2y = -p1.y + (p2.y - p0.y) * 0.1;
                                pathData += ` C ${cp1x} ${cp1y} ${cp2x} ${cp2y} ${p1.x} ${-p1.y}`;
                              } else {
                                // 중간 곡선: 연속성 보장
                                const cp1x = p0.x + (p1.x - traces[i-2].x) * 0.15;
                                const cp1y = -p0.y - (p1.y - traces[i-2].y) * 0.15;
                                const cp2x = p1.x - (p2.x - p0.x) * 0.15;
                                const cp2y = -p1.y + (p2.y - p0.y) * 0.15;
                                pathData += ` C ${cp1x} ${cp1y} ${cp2x} ${cp2y} ${p1.x} ${-p1.y}`;
                              }
                            }
                          }
                          
                          const color = radiusColors[padRadius];
                          const opacity = 0.7 - (rotIndex * 0.1);
                          
                          return (
                            <g key={`radius-${padRadius}-rotation-${rotIndex}`}>
                              <path 
                                d={pathData} 
                                fill="none" 
                                stroke={color} 
                                strokeWidth="2.5" 
                                opacity={Math.max(opacity, 0.3)}
                              />
                              {traces.length > 0 && (
                                <g>
                                  <circle 
                                    cx={traces[traces.length - 1].x} 
                                    cy={-traces[traces.length - 1].y} 
                                    r="4" 
                                    fill={color} 
                                    stroke="white" 
                                    strokeWidth="1" 
                                  />
                                  {traces.filter((_, i) => i % 3 === 0).slice(-8).map((point, i) => (
                                    <circle 
                                      key={i} 
                                      cx={point.x} 
                                      cy={-point.y} 
                                      r="1.5" 
                                      fill={color} 
                                      opacity={0.4 + i * 0.1} 
                                    />
                                  ))}
                                </g>
                              )}
                            </g>
                          );
                        });
                      })}
                      
                      <g transform="translate(-400, -350)">
                        <rect x="0" y="0" width="140" height={selectedPadRadii.length * 15 + 20} fill="white" opacity="0.95" stroke="gray" strokeWidth="1" rx="3"/>
                        <text x="5" y="15" fontSize="10" fontWeight="bold" fill="black">패드반지름별 자취</text>
                        {selectedPadRadii.map((radius, i) => {
                          const radiusColors = {
                            50: '#FF6B6B', 100: '#4ECDC4', 150: '#45B7D1', 200: '#FFA726',
                            250: '#66BB6A', 300: '#AB47BC', 350: '#FF7043', 390: '#8D6E63'
                          };
                          return (
                            <g key={radius}>
                              <line x1="5" y1={30 + i * 15} x2="15" y2={30 + i * 15} stroke={radiusColors[radius]} strokeWidth="3"/>
                              <text x="20" y={35 + i * 15} fontSize="8" fill="black">{radius}mm</text>
                            </g>
                          );
                        })}
                      </g>
                    </svg>
                  </div>
                </div>
              )}
            </div>
          </div>
        )}


      </div>
      
      {results && (
        <div className="bg-green-50 border border-green-200 rounded-lg p-4">
          <h3 className="text-lg font-semibold text-green-800 mb-3">시뮬레이션 완료</h3>
          <div className="grid grid-cols-1 sm:grid-cols-3 gap-4">
            <div className="bg-white p-3 rounded">
              <h4 className="font-medium text-green-700 mb-2">통계</h4>
              <ul className="text-sm space-y-1">
                <li>• 최대: {results.statistics.maxDistance.toFixed(1)} mm</li>
                <li>• 평균: {results.statistics.avgDistance.toFixed(1)} mm</li>
                <li>• 점 개수: {results.wafer_points.length}개</li>
              </ul>
            </div>
            <div className="bg-white p-3 rounded">
              <h4 className="font-medium text-green-700 mb-2">활성 웨이퍼</h4>
              <ul className="text-sm space-y-1">
                <li>• L웨이퍼: {activeWafers.left ? '✅' : '❌'}</li>
                <li>• R웨이퍼: {activeWafers.right ? '✅' : '❌'}</li>
              </ul>
            </div>
            <div className="bg-white p-3 rounded">
              <h4 className="font-medium text-green-700 mb-2">표시 그래프</h4>
              <ul className="text-sm space-y-1">
                <li>• 패드좌표계: {visibleGraphs.padCoordinate ? '✅' : '❌'}</li>
                <li>• 히트맵: {visibleGraphs.heatmap ? '✅' : '❌'}</li>
                <li>• 속도벡터: {visibleGraphs.velocityVector ? '✅' : '❌'}</li>
                <li>• 패드자취: {visibleGraphs.pad200mmTrace ? '✅' : '❌'}</li>
              </ul>
            </div>
          </div>
        </div>
      )}
    </div>
  );
};

export default WaferSimulator;